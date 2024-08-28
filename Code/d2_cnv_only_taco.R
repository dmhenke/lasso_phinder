# Load R libs ####
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


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(X_cnv)


# Run for list of genes on D2 ####
run_analysis <- function(y, gene){
  
  y <- y[!is.na(y)]
  
  # Find overlapping cell lines ####
  ok_cells <- intersect(names(y), rownames(X_cnv))
  
  X_cnv_ok  <- X_cnv[ok_cells, ]
  X_cnv_ok <- X_cnv_ok[, apply(X_cnv_ok, 2, var) > 0]
  
  y <- y[ok_cells]
  
  
  # Run LASSO ####
  scores <- get_scores(gene, ppi)
  
  scale <- function(matr){
    apply(matr, 2, function(x)
      (x - mean(x))/sd(x))
  }
  
 results_cnv <- run_reg_lasso(
    scale(X_cnv_ok), y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  
 return(results_cnv)
}

mclapply(
  mc.cores = 32, d2_genes, function(gene){
    y <- demeter2[, gene]
    res <- try(run_analysis(y, gene))
    if(class(res) == "try-error") return(NULL)
    save(
      res,
      file = paste0(
        "/storage/thinc/git_repos/lasso_phinder/Outputs/cnv_d2/",
        gene, ".RData"))
})

