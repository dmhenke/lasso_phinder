# Load R libs ####
library(caret)
library(data.table)
library(ggplot2)
library(glmnet)
library(parallel)


# Set working directory on taco ####
if(Sys.info()[[6]]=="Dafydd"){ setwd('../Code/')}else setwd("/storage/thinc/git_repos/lasso_phinder/Code")


# Source code ####
source("allfunctions.R")


# Load STRING data ####
load("../Data/ppi_w_symbols.RData")


# Load DepMap data ####
load("../Data/global.RData")


# Read skewness genes ####
d2_genes <- scan(
  "../Data/skewed_d2_genes.txt",
  what = character(), sep = "\n")

chr_genes <- scan(
  "../Data/skewed_chr_genes.txt",
  what = character(), sep = "\n")


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


# Run for list of genes on D2 ####
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
  
  scale <- function(matr){
    apply(matr, 2, function(x)
      (x - mean(x))/sd(x))
  }
  
  results_rna <- run_reg_lasso(
    scale(X_rna_ok), y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  results_cnv <- run_reg_lasso(
    scale(X_cnv_ok), y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  results_mut <- run_reg_lasso(
    scale(X_mut_ok), y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  # Combine non-zero coefficients into one X ####
  betas <- list()
  if(length(results_cnv) == 3){
    betas$CNV <- X_cnv_ok[, results_cnv$betas$betas_pen != 0]
    if(class(betas$CNV)[1] == "numeric"){
      betas$CNV <- matrix(betas$CNV, ncol = 1)
      colnames(betas$CNV) <- rownames(results_cnv$betas)[results_cnv$betas$betas_pen != 0]
    }
  }
  if(length(results_rna) == 3){
    betas$RNA <- X_rna_ok[, results_rna$betas$betas_pen != 0]
    if(class(betas$RNA)[1] == "numeric"){
      betas$RNA <- matrix(betas$RNA, ncol = 1)
      colnames(betas$RNA) <- rownames(results_rna$betas)[results_rna$betas$betas_pen != 0]
    }
  } 
  if(length(results_mut) == 3){
    betas$Mut <- X_mut_ok[, results_mut$betas$betas_pen != 0]
    if(class(betas$Mut)[1] == "numeric"){
      betas$Mut <- matrix(betas$Mut, ncol = 1)
      colnames(betas$Mut) <- rownames(results_mut$betas)[results_mut$betas$betas_pen != 0]
    }
  } 
  
  X_combined <- do.call(cbind, betas)
  
  genes <- unlist(lapply(betas, function(x)
    colnames(x)))
  omic <- unlist(lapply(names(betas), function(x){
    rep(x, ncol(betas[[x]])) 
  }))
  
  
  # Run LASSO again ####
  X <- scale(X_combined)
  
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


resl <- mclapply(mc.cores=32,c("demeter2","kronos"),function(scr){
#resl <- lapply(c("demeter2","kronos"),function(scr){
  if(scr=='kronos'){
    genes_test <- scan(
      "../Data/skewed_chr_genes.txt",
      what = character(), sep = "\n")
  }else if(scr=='demeter2'){
    genes_test <- scan(
      "../Data/skewed_d2_genes.txt",
      what = character(), sep = "\n")
  }
  # loop
  lapply(genes_test, function(x,whichY=scr){
    print(paste("Working on", x))
    if(whichY=="demeter2")y <- demeter2[, x] else y <- kronos[, x]
    y <- y[!is.na(y)]
    # check if analysis already run
    output_file <- paste0("../Outputs/",whichY,"_",x,"_results_multiomic.RData")
    if(!basename(output_file) %in% list.files("../Outputs/")){
      res <- try(run_analysis(y, x))
      #write results
	if(class(res)!="try-error"){
		save(res, file = output_file)		
	}
    }
    
  })
})
