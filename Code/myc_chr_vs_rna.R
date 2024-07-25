# Load R libs ####
library(caret)
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


# Define outcome ####
gene <- "MYC"
y <- kronos[, gene]
y <- y[!is.na(y)]


# Run LASSO ####
ok_cells <- intersect(
  rownames(X),
  names(y)
)

X <- X[ok_cells, ]
y <- y[ok_cells]

scores <- get_scores(gene, ppi)

results <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))


# Show how phi was inferred ####
tmp <- results$correlations
asplit <- split(1:nrow(tmp), tmp$run)
tmp <- do.call(cbind, lapply(asplit, function(x) tmp$cor[x]))
tmp <- apply(tmp, 2, function(x)
  (x - mean(x))/sd(x))

phi_range <- results$correlations$phi[1:30]
median_cor <- apply(tmp, 1, mean)

aframe <- data.frame(
  phi = phi_range, 
  cor = median_cor)

ggplot(aframe, aes(phi, cor)) +
  labs(
    title = "MYC [Chronos] vs RNA expression",
    y = "Average standardized correlation",
    x = "phi"
  ) +
  geom_line() + geom_point() +
  geom_vline(xintercept = results[[1]]$best_cor_phi) +
  theme_classic()


# Show non-zero coefficients ####
aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  cor = cor(X, y)[,1])

aframe$gene <- factor(
  aframe$gene, 
  levels = aframe$gene[order(aframe$betas_pen)])

ggplot(aframe[aframe$betas_pen != 0, ], aes(betas_pen, gene)) +
  labs(
    y = "Informative & relevant features",
    x = "Regularized LASSO coefficient"
  ) +
  geom_bar(stat = "identity") +
  theme_classic()


# Compare to regular correlation coefficients ####
aframe$label <- aframe$gene
aframe$label[aframe$betas_pen == 0] <- NA

ggplot(aframe, aes(cor, betas_pen,
                   label = label)) +
  labs(
    y = "Regularized LASSO coefficient",
    x = "Correlation coefficient"
  ) +
  geom_point() +
  ggrepel::geom_label_repel(max.overlaps = 20) +
  theme_classic()

plot_gene <- function(gene, y){
  subm <- data.frame(
    gene = X[, gene],
    y
  )  
  
  ggplot(subm, aes(gene, y)) +
    labs(
      x = paste(gene, "RNA levels"),
      y = "MYC dependency [Chronos]"
    ) +
    geom_smooth(method = "lm") +
    geom_point() +
    ggpubr::stat_cor() +
    theme_classic() 
}
plot_gene("MYC", y)
plot_gene("UBR5", y)
plot_gene("WEE1", y)

# Run one more time ####
results_new <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

aframe <- data.frame(
  results$betas,
  results_new$betas
)
ggplot(aframe,
       aes(betas_pen, betas_pen.1)) +
  labs(
    title = "Regularized betas",
    x = "Run 1",
    y = "Run 2"
  ) +
  geom_point() +
  ggpubr::stat_cor() +
  theme_classic()


# Run LASSO with MYC score set to 0 ####
scores["MYC"] <- 0

results_myc0 <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  results_myc0$betas
)

aframe$label <- aframe$gene
aframe$label[aframe$betas_pen == 0] <- NA

ggplot(aframe,
       aes(betas_pen, betas_pen.1, label = label)) +
  labs(
    title = "Regularized betas",
    x = "Original",
    y = "MYC score set to 0"
  ) +
  geom_point() +
  ggrepel::geom_label_repel() +
  ggpubr::stat_cor() +
  theme_classic()


# Compare biomarkers to predict MYC dependency in D2 ####
myc_d2 <- demeter2[, "MYC"]
myc_d2 <- myc_d2[!is.na(myc_d2)]

X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))

ok_cells <- intersect(
  rownames(X),
  names(myc_d2)
)

X <- X[ok_cells, ]
myc_d2 <- myc_d2[ok_cells]

genes <- aframe$gene[aframe$betas_pen != 0]
genes_cor <- tail(aframe$gene[order(abs(aframe$cor))], length(genes))

library(ggVennDiagram)
ggVennDiagram(list(genes, genes_cor))
intersect(genes, genes_cor)

fitControl <- trainControl(
  method = "cv", number = 10, savePredictions = T)

subm <- data.frame(myc_d2, X[, setdiff(genes, genes_cor)])
afit_lasso <- train(
  myc_d2 ~ ., data = subm, method = "lm", trControl = fitControl)

subm <- data.frame(myc_d2, X[, setdiff(genes_cor, genes)])
afit_cor <- train(
  myc_d2 ~ ., data = subm, method = "lm", trControl = fitControl)

aframe <- data.frame(
  lasso = afit_lasso$resample$Rsquared,
  cor = afit_cor$resample$Rsquared)
aframe <- reshape2::melt(aframe)

ggplot(aframe, aes(variable, value)) +
  labs(
    title = "Predicting MYC dependency in D2",
    y = "R-squared",
    x = "Biomarkers from"
  ) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.1) +
  ggpubr::stat_compare_means(method = "t.test") +
  theme_classic()

