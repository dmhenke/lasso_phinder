# Load R libs ####
library(caret)
library(data.table)
library(ggplot2)
library(glmnet)
library(parallel)


# Set working directory on taco ####
if(Sys.info()[[6]]=="Dafydd"){ setwd("../");setwd(paste0(getwd(),'/Code/'))}else setwd("/storage/thinc/git_repos/lasso_phinder/Code")


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
# X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]
# offsetlog2 <- max(log2(cnv[which(cnv>2)]))+0.2
# X_cnv[which(cnv<2)] <- cnv[which(cnv<2)] +2^(-offsetlog2)*(2-cnv[which(cnv<2)] )


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut
X_mut <- X_mut[, colSums(X_mut) >= 5]


# resl <- mclapply(mc.cores=32,c("demeter2","kronos"),function(scr){
resl <- lapply(c("demeter2","kronos"),function(scr){
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
      res <- run_analysis(y, x)
      #write results
      save(res, file = output_file)
    }
  })
})


#single omic
# EGFR & D2 for CNV
# res_sing <- lapply(c("demeter2","kronos")[1],function(scr,which_omic="CNV3"){
  res_sing <- lapply(c("demeter2"),function(scr,which_omic="CNVskew2"){
    genes_test <-"EGFR"
  #loop
  lapply(genes_test, function(x,whichY=scr){
    print(paste("Working on", x))
    if(whichY=="demeter2")y <- demeter2[, x] else y <- kronos[, x]
    y <- y[!is.na(y)]
    # check if analysis already run
    output_file <- paste0("../Outputs/",whichY,"_",x,"_results_",which_omic,"_omic.RData")
    if(!basename(output_file) %in% list.files("../Outputs/")){
      res <- run_analysisSingle(y, gene=x,omics=list(X_cnv),which_omic=c("CNV"))

      #write results
      save(res, file = output_file)
    }
  })
})

# PRMT5
  system.time({
  run_analysisSingle(whichY="demeter2", gene='E',omics=list(X_cnv),which_omic="CNV")
  })
  # "Lambda min: 0.0212"
  # "Best phi based on RMSE: 0.069"
  
  # Oct4 (Pou5f1), Sox2, Klf4 and cMyc.
  
  lapply(c('POU5F1','SOX2','KLF4','MYC'),function(x)  run_analysisSingle(whichY="demeter2", gene=x,omics=list(X_cnv),which_omic="CNV"))
  
  