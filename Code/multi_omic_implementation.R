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


# Define function ####
run_analysis <- function(y, gene){
  
  # Find overlapping cell lines ####
  ok_cells <- intersect(names(y), rownames(X_rna))
  ok_cells <- intersect(ok_cells, rownames(X_cnv))
  ok_cells <- intersect(ok_cells, rownames(X_mut))
  
  
  # Remove features without variance ####
  X_rna_ok  <- X_rna[ok_cells, ]
  X_rna_ok <- X_rna_ok[, apply(X_rna_ok, 2, var) > 0]
  
  X_mut_ok <- X_mut[ok_cells, ]
  X_mut_ok <- X_mut_ok[, colSums(X_mut_ok) >= 5]
  
  X_cnv_ok  <- X_cnv[ok_cells, ]
  X_cnv_ok <- X_cnv_ok[, apply(X_cnv_ok, 2, var) > 0]
  
  y <- y[ok_cells]
  
  
  # Run LASSO ####
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
  
  # Combine non-zero coefficients into one X ####
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
  
  
  # Run LASSO again ####
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
  
  aframe <- data.frame(
    gene = genes,
    omic,
    correl = cor(X, y),
    betas,
    betas_pen,
    score = scores[genes])
  
  return(aframe)
}

gene <- "VPS13D"
y <- kronos[, gene]
out <- run_analysis(y, gene)

# ####
res <- res_d2
bad <- which(unlist(lapply(res, class)) == "try-error")
good <- names(res)[-bad]
d2_combined <- do.call(rbind, lapply(good, function(x)
  data.frame(target = x, res[[x]])))

tmp <- d2_combined
tmp <- tmp[tmp]


res <- res_chr
bad <- which(unlist(lapply(res, class)) == "try-error")
good <- names(res)[-bad]
d2_combined <- do.call(rbind, lapply(good, function(x)
  data.frame(target = x, res[[x]])))

