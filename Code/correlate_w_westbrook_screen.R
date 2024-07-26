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


# Run biolasso ####
gene <- "MYC"
y <- kronos[, gene]

ok_cells <- intersect(names(y), rownames(X_rna))
ok_cells <- intersect(ok_cells, rownames(X_cnv))
ok_cells <- intersect(ok_cells, rownames(X_mut))

X_rna_ok  <- X_rna[ok_cells, ]
X_rna_ok <- X_rna_ok[, apply(X_rna_ok, 2, var) > 0]

X_mut_ok <- X_mut[ok_cells, ]
X_mut_ok <- X_mut_ok[, colSums(X_mut_ok) >= 5]

X_cnv_ok  <- X_cnv[ok_cells, ]
X_cnv_ok <- X_cnv_ok[, apply(X_cnv_ok, 2, var) > 0]

y <- y[ok_cells]
  
scores <- get_scores(gene, ppi)

results_rna <- run_reg_lasso(
  X_rna_ok, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_cnv <- run_reg_lasso(
  X_cnv_ok, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_mut <- run_reg_lasso(
  X_mut_ok, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))
  
betas <- list()
if(length(results_cnv) == 3){
  betas$CNV <- X_cnv_ok[, results_cnv$betas$betas_pen != 0]
  if(class(betas$CNV)[1] == "numeric"){
    betas$CNV <- matrix(betas$CNV, ncol = 1)
  }
}
if(length(results_rna) == 3){
  betas$RNA <- X_rna_ok[, results_rna$betas$betas_pen != 0]
  if(class(betas$RNA)[1] == "numeric"){
    betas$RNA <- matrix(betas$RNA, ncol = 1)
  }
} 
if(length(results_mut) == 3){
  betas$Mut <- X_mut_ok[, results_mut$betas$betas_pen != 0]
  if(class(betas$Mut)[1] == "numeric"){
    betas$Mut <- matrix(betas$Mut, ncol = 1)
  }
} 
  
X_combined <- do.call(cbind, betas)

genes <- unlist(lapply(betas, colnames))
omic <- unlist(lapply(names(betas), function(x){
  rep(x, ncol(betas[[x]])) 
}))

X <- X_combined

lambda_min <- find_lambda(X, y, plot = F)

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

out <- data.frame(
  gene = genes,
  omic,
  correl = cor(X, y),
  betas,
  betas_pen,
  score = scores[genes])


# Load screen results ####
screen <- readxl::read_excel(
  "/Users/lukas/Downloads/2020.12.22.SYNTHETIC.LETHALS.gene_table_integrated_v2.TW3.xlsx", 
  sheet = "tam_unt")
screen <- screen[, c("gene", "pool", "ptpn12_tam_unt_l2fc", "ptpn12_tam_unt_pval")]
screen <- screen[sort.list(screen$ptpn12_tam_unt_pval), ]

genes <- unique(out$gene[out$betas_pen != 0])
genes <- intersect(genes, screen$gene)

screen$lasso <- screen$gene %in% genes
ggplot(
  screen[order(screen$lasso), ],
  aes(ptpn12_tam_unt_l2fc, -log10(ptpn12_tam_unt_pval),
      color = lasso)) +
  labs(
    title = "LASSO results vs. screen",
    x = "Fold change",
    y = "-log10 P-value") +
  geom_point() +
  ggrepel::geom_text_repel(
    data = screen[screen$lasso, ][1:20, ],
    aes(label = gene)) +
  theme_classic()


genes <- rownames(results_cnv$betas)[results_cnv$betas$betas_pen != 0]
genes <- intersect(genes, screen$gene)

subm <- data.frame(
  screen[match(genes, screen$gene),],
  beta = results_cnv$betas[match(genes, rownames(results_cnv$betas)), ]
)

ggplot(
  subm,
  aes(ptpn12_tam_unt_l2fc, -log10(ptpn12_tam_unt_pval),
      color = beta.betas_pen > 0,
      label = gene)) +
  labs(
    title = "LASSO results vs. screen",
    x = "Fold change",
    y = "-log10 P-value") +
  geom_point() +
  ggrepel::geom_text_repel() +
  theme_classic()

fisher.test(table(
  subm$ptpn12_tam_unt_l2fc < 0, subm$beta.betas_pen > 0))

