# Load R libs ####
library(ggplot2)


# Source code ####
source("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Code/allfunctions.R")


# Load data ####
load("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/ppi_w_symbols.RData")
load("/Users/lukas/OneDrive/Miko/THINC/projects/TCGA/data_for_lasso.RData")


# Calculate gene-wise CNA burden ####
meta$total_burden <- apply(gistic_merged, 2, function(x) mean(x != 0))


# Test some associations ####
aframe <- data.frame(
  meta,
  cnv  = gistic_merged["U2AF2", ] < 0)

coefficients(summary(lm(
  u2af2_csj ~ Project.ID + gii + gender + cnv,
  data = aframe)))

coefficients(summary(lm(
  u2af2_csj ~ Project.ID + total_burden + gender + cnv,
  data = aframe)))


# Restrict to BRCA ####
brca <- which(meta$Project.ID == "TCGA-BRCA" &
                !is.na(meta$gii) & !is.na(meta$dhx15_csj))
cnv_brca <- gistic_merged[, brca]
chr14_brca <- cnv_brca#[which(gene_list$updated == 14), ]
meta_brca <- meta[brca, ]

#resids <- residuals(lm(dhx15_csj ~ gii + gender, data = meta_brca))
resids <- residuals(lm(u2af2_csj ~ gii + gender, data = meta_brca))
y <- resids

gene <- "U2AF2"
scores <- get_scores(gene, ppi)

X <- t(chr14_brca)
X <- apply(X, 2, function(x) as.numeric(x < 0))

results <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))


# Create plots ####
tmp <- data.frame(
  gene = colnames(X), 
  results$betas,
  correl = cor(X, y)[,1],
  chr = gene_list$updated
)

tmp$rank <- 1:nrow(tmp)

tmp$label <- colnames(X)
tmp$label[which(abs(tmp$betas_pen) < 1e-5)] <- NA

ggplot(tmp, aes(x = rank, y = correl,
                label = label, 
                size = abs(betas_pen), 
                color = (gene == "U2AF2"))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.5) +
  ggrepel::geom_text_repel(max.overlaps = 50, color = "black") +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Rank", y = "Correlation with residuals")

