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
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(log2(X_cnv))


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut
X_mut <- X_mut[, colSums(X_mut) >= 5]


# Read skewness genes ####
d2_genes <- scan(
  "../Data/skewed_d2_genes.txt",
  what = character(), sep = "\n")

chr_genes <- scan(
  "../Data/skewed_chr_genes.txt",
  what = character(), sep = "\n")


# Load LASSO results ####
load("../Outputs/d2_results_multiomic.RData")
res_d2 <- res
names(res_d2) <- d2_genes

load("../Outputs/chr_results_multiomic.RData")
res_chr <- res
names(res_chr) <- chr_genes



# ####
res <- res_d2
bad <- which(unlist(lapply(res, class)) == "try-error")
good <- names(res)[-bad]
d2_combined <- do.call(rbind, lapply(good, function(x)
  data.frame(target = x, res[[x]])))

tmp <- d2_combined
tmp <- tmp[tmp$betas_pen != 0, ]

subm <- tmp[tmp$target == "AURKA" & tmp$omic == "CNV", ]

ggplot(subm, aes(correl, betas_pen)) +
  geom_point() +
  labs(
    x = "Correlation coefficient",
    y = "Regularized LASSO beta"
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(aes(label = gene)) +
  theme_classic()

y <- demeter2[, "AURKA"]
y <- y[!is.na(y)]

aframe <- data.frame(
  cnv = X_cnv[match(names(y), rownames(X_cnv)), c("AURKA", "TACC2")],
  mut = X_mut[match(names(y), rownames(X_mut)), c("PASK", "SNX31")],
  y  
)
ggplot(aframe, aes(cnv.AURKA, y)) +
  geom_point() +
  stat_cor() +
  theme_classic()


res <- res_chr
bad <- which(unlist(lapply(res, class)) == "try-error")
good <- names(res)[-bad]
d2_combined <- do.call(rbind, lapply(good, function(x)
  data.frame(target = x, res[[x]])))

