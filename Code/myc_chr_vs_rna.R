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
median_cor <- apply(tmp, 1, mean)

aframe <- data.frame(
  phi = seq(0, 1, length = 30), 
  cor = median_cor)

ggplot(aframe, aes(phi, cor)) +
  labs(
    title = "MYC [Chronos] vs RNA expression",
    y = "Average z-scored RMSE",
    x = "phi"
  ) +
  geom_line() + geom_point() +
  geom_vline(xintercept = results[[1]]) +
  theme_classic()

ggsave(
  "../Outputs/figs/myc_rna_phi_curve.pdf",
  width = 5, height = 5)


# Show non-zero coefficients ####
aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  cor = cor(X, y)[,1])

aframe$gene <- factor(
  aframe$gene, 
  levels = aframe$gene[order(aframe$betas_pen)])

aframe <- aframe[order(-abs(aframe$betas_pen)), ]

ggplot(aframe[aframe$betas_pen != 0, ][1:40, ], aes(betas_pen, gene)) +
  labs(
    y = "Informative & relevant features [top 40]",
    x = "Regularized LASSO coefficient"
  ) +
  geom_bar(stat = "identity") +
  theme_classic()

ggsave(
  "../Outputs/figs/myc_rna_top40_barplot.pdf",
  width = 5, height = 8)


# Compare to standard LASSO coefficients ####
aframe$label <- aframe$gene
aframe$label[aframe$betas_pen == 0] <- NA

p1 <- ggplot(aframe, aes(betas, betas_pen,
                   label = label)) +
  labs(
    y = "Regularized LASSO coefficient",
    x = "LASSO coefficient"
  ) +
  geom_point() +
  ggrepel::geom_label_repel(max.overlaps = 20) +
  theme_classic()

p2 <- ggplot(aframe[aframe$gene != "MYC", ], aes(betas, betas_pen, 
                   color = cor,
                   label = label)) +
  labs(
    y = "Regularized LASSO coefficient",
    x = "LASSO coefficient"
  ) +
  scale_color_gradient2(
    low = "blue", mid = "grey", high = "red",
    breaks = seq(-0.2, 0.2, length = 5)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  theme_classic()

p <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(
  plot = p,
  filename = "../Outputs/figs/myc_betas_reg_vs_standard.pdf",
  width = 11, height = 5)


# Highlight some genes ####
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
p1 <- plot_gene("MYC", y)
p2 <- plot_gene("STAT5A", y)
p3 <- plot_gene("SMARCB1", y)

p <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3)

ggsave(
  plot = p, 
  filename = "../Outputs/figs/myc_single_gene_examples.pdf",
  width = 15, height = 5)


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

ggsave(
  "../Outputs/figs/myc_reproduce.pdf",
  width = 5, height = 5)


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

ggsave(
  "../Outputs/figs/myc_score_set_to_zero.pdf",
  width = 5, height = 5)


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


