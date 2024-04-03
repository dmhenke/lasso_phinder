# Load R libs ####
library(caret)
library(data.table)
library(ggplot2)
library(glmnet)


# Source code ####
source("allfunctions.R")


# Load STRING data ####
load("../Data/ppi_w_symbols.RData")


# Load DepMap data ####
load("../Data/global.RData")


# Normalize RNA expression data ####
X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(X_cnv)


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut


# Define outcome ####
gene <- "PKMYT1"
y <- demeter2[, gene]
y <- kronos[, gene]
y <- y[!is.na(y)]

tissue <- sample_info[names(y), "lineage"]
y_resid <- residuals(lm(y ~ tissue))
y <- y_resid


# Find overlapping cell lines ####
ok_cells <- intersect(names(y), rownames(X_rna))
ok_cells <- intersect(ok_cells, rownames(X_cnv))
ok_cells <- intersect(ok_cells, rownames(X_mut))


# Remove features without variance ####
X_rna  <- X_rna[ok_cells, ]
X_rna <- X_rna[, apply(X_rna, 2, var) > 0]

X_mut <- X_mut[ok_cells, ]
X_mut <- X_mut[, colSums(X_mut) >= 5]

X_cnv  <- X_cnv[ok_cells, ]
X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]

y <- y[ok_cells]


# Run LASSO ####
scores <- get_scores(gene, ppi)

results_rna <- run_reg_lasso(
  X_rna, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_cnv <- run_reg_lasso(
  X_cnv, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_mut <- run_reg_lasso(
  X_mut, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))


# Combine non-zero coefficients into one X ####
X_combined <- cbind(
  X_rna[, results_rna$betas$betas_pen != 0],
  X_cnv[, results_cnv$betas$betas_pen != 0],
  X_mut[, results_mut$betas$betas_pen != 0])

genes <- c(rownames(results_rna$betas)[results_rna$betas$betas_pen != 0],
           rownames(results_cnv$betas)[results_cnv$betas$betas_pen != 0],
           rownames(results_mut$betas)[results_mut$betas$betas_pen != 0])

omic <- c(rep("rna", sum(results_rna$betas$betas_pen != 0)),
          rep("cnv", sum(results_cnv$betas$betas_pen != 0)),
          rep("mut", sum(results_mut$betas$betas_pen != 0)))

X <- X_combined

lambda_min <- find_lambda(X, y, plot = T)

afit <- glmnet(
  X, y, 
  alpha = 1, 
  lambda = lambda_min)
betas <- afit$beta[,1]

penalties <- scores[match(genes, names(scores))]
names(penalties) <- genes
penalties[is.na(penalties)] <- 0
penalties <- penalties/max(scores)

correls <- get_rmse_from_xvalid(
  X, y,
  penalties,
  lambda_min = lambda_min,
  phi_range = seq(0, 1, length = 30), n_folds = 10)

best_phi <- find_best_phi_rmse(
  correls,
  phi_range = seq(0, 1, length = 30))

afit <- glmnet(
  X,
  y,
  alpha = 1,
  lambda = lambda_min,
  penalty.factor = 1 - penalties * best_phi)
betas_pen <- afit$beta[,1]

aframe <- data.frame(
  gene = genes,
  omic, 
  betas,
  betas_pen, score = scores[genes])

ggplot(aframe, aes(x = betas, 
                   y = betas_pen, 
                   color = omic,
                   label = gene)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  theme_classic()

