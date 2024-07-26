# Load R libs ####
library(GEOquery)
library(ggplot2)
library(preprocessCore)


# Load data ####
load("/storage/thinc/projects/westbrook_lab/dhx15_paper/analysis_dhx15_manuscript/outputs/global.RData")

drug <- read.csv("/storage/thinc/projects/resources/depmap/data/drug_sensitivity/depmap_drug_auc.csv", row.names = 1)
drug_meta <- read.csv("/storage/thinc/projects/resources/depmap/data/drug_sensitivity/depmap_drug_metadata.csv", row.names = 1)


# Normalize RNA expression data ####
X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(log2(X_cnv))


# Run biolasso ####
source("/storage/thinc/git_repos/lasso_phinder/Code/allfunctions.R")
load("/storage/thinc/git_repos/lasso_phinder/Data/ppi_w_symbols.RData")

scores <- get_scores("EGFR", ppi)

drug_meta[drug_meta$drug_name == "ERLOTINIB", ]

y <- drug[, 1086]
names(y) <- rownames(drug)

y <- demeter2[, "EGFR"]
y <- y[!is.na(y)]

ok <- intersect(rownames(X_rna), names(y))

results_rna <- run_reg_lasso(
  scale(X_rna[ok,]), y[ok], scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_rna <- data.frame(
  gene = colnames(X_rna),
  results_rna$betas,
  cor = cor(X_rna[ok,], y[ok])[,1],
  score = scores[match(colnames(X_rna), names(scores))]
)


plot(X_rna[ok, "ERBB3"], y[ok])
plot(X_rna[ok, "MAPK14"], y[ok])

ok <- intersect(rownames(X_cnv), names(y))

results_cnv <- run_reg_lasso(
  scale(X_cnv[ok,]), y[ok], scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

results_cnv <- data.frame(
  gene = colnames(X_cnv),
  results_cnv$betas,
  cor = cor(X_cnv[ok,], y[ok])[,1],
  score = scores[match(colnames(X_cnv), names(scores))]
)


# Get public available data ####
gset <- getGEO("GSE31625", GSEMatrix =TRUE, getGPL=TRUE)[[1]]
ex <- exprs(gset)
ex <- log2(ex)

norm <- normalize.quantiles(ex)

treat <- gset$`erlotinib sensitivity:ch1`
gene_info <- gset@featureData@data
genes <- gene_info$`Gene Symbol`

gene_info[which(genes == "ERBB3"),]


plot_gene <- function(gene){
  if(length(which(genes == gene)) == 1){
    subm <- data.frame(
      treat,
      norm[which(genes == gene),])
  }
  if(length(which(genes == gene)) > 1){
    subm <- data.frame(
      treat,
      t(norm[which(genes == gene),]))
  }
  
  subm <- reshape2::melt(subm, id.vars = "treat")
  ggplot(subm, aes(treat, value)) +
    facet_wrap(~ variable, scales = "free_y", nrow = 1) +
    geom_boxplot() +
    theme_classic()  
}
plot_gene("EGF")

res <- t(apply(norm, 1, function(x)
  coefficients(summary(lm(x ~ treat)))[2, c(1, 4)]))
res <- data.frame(res)
colnames(res) <- c("coef", "pval")
res$gene <- gene_info$`Gene Symbol`

res <- res[sort.list(res$pval), ]

res$label <- NA
res$label[1:50] <- res$label[1:50]

down <- head(results_rna$gene[order(results_rna$cor)], 500)
up <- tail(results_rna$gene[order(results_rna$cor)], 500)

tmp <- split(results_rna$gene[results_rna$betas != 0],
              results_rna$cor[results_rna$betas != 0] > 0)
tmp <- split(results_rna$gene[results_cnv$betas_pen != 0],
             results_rna$cor[results_cnv$betas_pen != 0] > 0)
tmp <- split(results_cnv$gene[results_cnv$betas_pen != 0],
             results_cnv$cor[results_cnv$betas_pen != 0] > 0)
down <- tmp[[1]]
up <- tmp[[2]]
#up <- results_rna$gene[results_rna$betas_pen > 0]

res$class <- "none"
res$class[res$gene %in% down] <- "down"
res$class[res$gene %in% up] <- "up"

ggplot(res[res$class != "none", ],
       aes(coef, -log10(pval), label = gene, color = class)) +
  geom_point() +
  geom_text() +
  theme_classic()

