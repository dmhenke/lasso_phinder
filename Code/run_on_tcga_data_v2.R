# Load R libs ####
library(ggplot2)


# Source code ####
source("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Code/allfunctions.R")


# Load data ####
load("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/ppi_w_symbols.RData")
load("/Users/lukas/OneDrive/Miko/THINC/projects/TCGA/data_for_lasso.RData")


pan <- read.csv(
  "/Users/lukas/Downloads/gene_loss_pancancer.csv")

ggplot(pan, aes(rank, association_score)) +
  geom_point() +
  theme_classic()


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

brca <- which(!is.na(meta$gii) & !is.na(meta$dhx15_csj) & !is.na(meta$gender))
brca <- sample(brca, 100)

cnv_brca <- gistic_merged[, brca]
chr14_brca <- cnv_brca#[which(gene_list$updated == 14), ]
meta_brca <- meta[brca, ]

#resids <- residuals(lm(u2af2_csj ~ gii + gender, data = meta_brca))
#gene <- "U2AF2"

gene <- "DHX15"
resids <- residuals(lm(dhx15_csj ~ Project.ID + gii + gender, data = meta_brca))

y <- resids

scores <- get_scores(gene, ppi)

X <- t(chr14_brca)
X <- apply(X, 2, function(x) as.numeric(x < 0))

system.time(results <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30)))


# Create plots ####
tmp <- data.frame(
  gene = colnames(X), 
  results$betas,
  score = scores[match(colnames(X), names(scores))],
  correl = cor(X, y)[,1]
)
aframe <- data.frame(
  pan,
  tmp[match(pan$gene, tmp$gene),])


ggplot(aframe, 
       aes(x = rank, y = association_score, color = factor(chr))) +
  labs(
    title = "TCGA",
    subtitle = "Pan-cancer",
    x = "Genes ordered by genomic location", 
    y = "Association score") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  scale_color_manual(
    values = rep(c("grey", "black"), 11)
  ) +
  theme_classic() +
  theme(legend.position = "none")


ggplot(aframe[aframe$chr == 3,],
       aes(x = rank, y = association_score, color = factor(chr))) +
      #aes(x = rank, y = correl, color = factor(chr))) +
  labs(
    title = "TCGA",
    subtitle = "Pan-cancer",
    x = "Genes ordered by genomic location (chr3)", 
    y = "Association score") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  scale_color_manual(
    values = rep(c("grey", "black"), 11)
  ) +
  theme_classic() +
  theme(legend.position = "none")


ggplot(aframe[aframe$chr == 3,],
       aes(x = rank, y = betas, color = factor(chr))) +
  labs(
    x = "Genes ordered by genomic location (chr3)", 
    y = "Baseline coefficient") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  scale_color_manual(
    values = rep(c("grey", "black"), 11)
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(
    data = aframe[aframe$chr == 3 &
                    aframe$betas != 0, ],
    aes(label = gene, size = abs(betas)),
    color = "red")


ggplot(aframe[aframe$chr == 3,],
       aes(x = rank, y = (-1)*betas_pen, color = factor(chr))) +
  labs(
    title = "TCGA",
    subtitle = "Pan-cancer",
    x = "Genes ordered by genomic location (chr3)", 
    y = "Bio-primed coefficient") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  scale_color_manual(
    values = rep(c("grey", "black"), 11)
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(
    data = aframe[aframe$chr == 3 &
                    aframe$betas_pen != 0, ],
    aes(label = gene, size = abs(betas_pen)),
    color = "red")

tmp <- aframe[aframe$chr == 3, ]
tmp <- reshape2::melt(tmp, measure.vars = c("betas", "betas_pen"))
ggplot(tmp,
       aes(x = rank, y = (-1)*value, color = factor(chr))) +
  facet_wrap(~ variable, nrow = 2) +
  labs(
    x = "Genes ordered by genomic location (chr3)", 
    y = "Bio-primed coefficient") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  scale_color_manual(
    values = rep(c("grey", "black"), 11)
  ) +
  theme_classic() +
  ylim(-1.6e-5, 1.6e-5) +
  theme(legend.position = "none")
