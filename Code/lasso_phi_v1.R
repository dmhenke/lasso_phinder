# Load R libs ####
library(biomaRt)
library(data.table)
library(ggplot2)
library(glmnet)


# Load String data ####
wrkfldr <- "/mnt/data/user/david/lasso/"

ppi <- fread("C:/Users/Dafydd/Documents/Projects/lasso_PPI/9606.protein.links.v12.0.txt")
ppi <- fread("/mnt/data/user/david/9606.protein.links.v12.0.txt")

ppi$protein1 <- gsub("9606.", "", ppi$protein1, fixed = T)
ppi$protein2 <- gsub("9606.", "", ppi$protein2, fixed = T)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol"), mart = mart)

ppi$gene1 <- gene_info$hgnc_symbol[match(ppi$protein1, gene_info$ensembl_peptide_id)]
ppi$gene2 <- gene_info$hgnc_symbol[match(ppi$protein2, gene_info$ensembl_peptide_id)]


# Load DepMap data ####
# load("C:/Users/Dafydd/Documents/Projects/lasso_PPI/global.RData")
load(paste0(wrkfldr,"/global.RData"))


# Define helper functions ####
get_scores <- function(gene, ppi){
  tmp <- ppi[(ppi$gene1 %in% c(gene) | 
                ppi$gene2 %in% c(gene)), ]
  tmp <- tmp[tmp$gene1 != "" & tmp$gene2 != "", ]
  scores <- as.numeric(tmp[["combined_score"]][tmp$gene1 == gene])
  names(scores) <- as.character(tmp[["gene2"]][tmp$gene1 == gene])
  scores[gene] <- 1000
  return(scores)
}

find_lambda <- function(X, y, plot = T){
  fitcv <- cv.glmnet(
    X, y, 
    alpha = 1, 
    nfolds = round(nrow(X)/30,0),
    lambda = NULL)
  if(plot) plot(fitcv, xvar = "lambda", label = T)
  # fitcv$lambda.min
  fitcv$lambda.1se # more conservative vs 'lambda.min'
}

find_best_phi <- function(correls, phi_range, plot = T){
  median_correls <- unlist(lapply(split(correls$cor,correls$phi), median))
  median_rmse <- unlist(lapply(split(correls$rmse,correls$phi), median))
  
  aframe <- data.frame(
    phi = phi_range, 
    cor = median_correls,
    rmse = median_rmse)
  
  # afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
  afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- predict(afit, aframe)
  diff_corr <- aframe$cor - preds
  afit2 <- lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds2 <- predict(afit2, aframe)
  diff_rmse <- aframe$rmse - preds2
  
  best_cor_phi <- phi_range[which(diff_corr == max(diff_corr))]
  best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]
  
  if(plot){
    # ggplot(aframe, aes(phi, cor)) +
    #   labs(
    #     x = 'phi',
    #     y = 'Correlation'
    #   ) +
    #   geom_vline(xintercept = best_phi, linetype = "dashed") +
    #   geom_point() +
    #   theme_classic()
    aframe$size=1;aframe$size[aframe$phi%in% c(best_cor_phi,best_rmse_phi)]<-2
    ggplot(aframe, aes(rmse, cor,color=phi,size=size)) +
      labs(
        x = 'RSME',
        y = 'Correlation'
      ) +
      geom_point() +
      geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
      theme_classic()
    diffame <-data.frame(diff_corr,diff_rmse,phi=phi_range,size=1)
    diffame$size=1;diffame$size[diffame$phi%in% c(best_cor_phi,best_rmse_phi)]<-2
    ggplot(diffame, aes(diff_rmse, diff_corr,color=phi,size=size)) +
      labs(
        x = 'diff RSME',
        y = 'diff Correlation'
      ) +
      geom_point() +
      geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
      theme_classic()
  }
  
  return(data.frame(best_cor_phi,best_rmse_phi))
}

run_reg_lasso <- function(X, y, scores,
                          n_folds = 10,
                          phi_range = seq(0, 1, length = 30)){
  lambda_min <- find_lambda(X, y, plot = F)
  
  # print(paste("Lambda min:", lambda_min))
  print(paste("Lambda 1st SE:", round(lambda_min,4)))
  
  afit <- glmnet(
    X, y, 
    alpha = 1, 
    lambda = lambda_min)
  betas <- afit$beta[,1]
  
  penalties <- scores[match(colnames(X), names(scores))]
  names(penalties) <- colnames(X)
  penalties[is.na(penalties)] <- 0
  penalties <- penalties/1000
  
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  
  correls <- do.call(rbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    
    do.call(rbind,lapply(phi_range, function(phi){
      lasso_tr <- glmnet(
        X[train,],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
      
      pred <- predict(lasso_tr,X[test,])
      rmse <- sqrt(apply((y[test]-pred)^2,2,mean))
      corout <- try(cor(pred, y[test]),silent = T)
      data.frame(cor=cor(pred, y[test]) ,rmse=rmse,run=x,phi=phi)
    }))
  }))
  
  best_phi <- find_best_phi(correls, phi_range)$best_cor_phi
  
  
  afit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi)
  betas_pen <- afit$beta[,1]
  
  best_phi <- find_best_phi(correls, phi_range, plot = F)$best_cor_phi
  
  print(paste("Phi best:", round(best_phi,4)))
  
  afit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi)
  betas_pen <- afit$beta[,1]
  
  return(data.frame(
    betas, 
    betas_pen))
}

run_reg_lasso_phionly <- function(X, y, scores,
                                  n_folds = 10,
                                  phi_range = seq(0, 1, length = 30)){
  lambda_min <- find_lambda(X, y, plot = F)
  
  # print(paste("Lambda min:", lambda_min))
  print(paste("Lambda 1st SE:", round(lambda_min,4)))
  
  afit <- glmnet(
    X, y, 
    alpha = 1, 
    lambda = lambda_min)
  betas <- afit$beta[,1]
  
  penalties <- scores[match(colnames(X), names(scores))]
  names(penalties) <- colnames(X)
  penalties[is.na(penalties)] <- 0
  penalties <- penalties/1000
  
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  
  correls <- do.call(rbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    
    do.call(rbind,lapply(phi_range, function(phi){
      lasso_tr <- glmnet(
        X[train,],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
      
      pred <- predict(lasso_tr,X[test,])
      rmse <- sqrt(apply((y[test]-pred)^2,2,mean))
      if(var(pred) == 0) cor_out <- NA else cor_out <- cor(pred, y[test])
      data.frame(cor= cor_out,rmse=rmse,run=x,phi=phi)
    }))
  }))
  if(any(is.na(correls))) data.frame(best_cor_phi=NA,best_rmse_phi=NA) else find_best_phi(correls, phi_range,plot=F)
}
wrap_regLassoPhionly <- function(gene){
  library(biomaRt)
  library(data.table)
  library(ggplot2)
  library(glmnet)
  
  get_scores <- function(gene, ppi){
    tmp <- ppi[(ppi$gene1 %in% c(gene) | 
                  ppi$gene2 %in% c(gene)), ]
    tmp <- tmp[tmp$gene1 != "" & tmp$gene2 != "", ]
    scores <- as.numeric(tmp[["combined_score"]][tmp$gene1 == gene])
    names(scores) <- as.character(tmp[["gene2"]][tmp$gene1 == gene])
    scores[gene] <- 1000
    return(scores)
  }
  
  find_lambda <- function(X, y, plot = T){
    fitcv <- cv.glmnet(
      X, y, 
      alpha = 1, 
      nfolds = round(nrow(X)/30,0),
      lambda = NULL)
    if(plot) plot(fitcv, xvar = "lambda", label = T)
    # fitcv$lambda.min
    fitcv$lambda.1se # more conservative vs 'lambda.min'
  }
  
  find_best_phi <- function(correls, phi_range, plot = T){
    median_correls <- unlist(lapply(split(correls$cor,correls$phi),median))
    median_rmse <- unlist(lapply(split(correls$rmse,correls$phi),median))
    
    aframe <- data.frame(
      phi = phi_range, 
      cor = median_correls,
      rmse = median_rmse)
    
    # afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
    afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
    preds <- predict(afit, aframe)
    diff_corr <- aframe$cor - preds
    afit2 <- lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
    preds2 <- predict(afit2, aframe)
    diff_rmse <- aframe$rmse - preds2
    
    best_cor_phi <- phi_range[which(diff_corr == max(diff_corr))]
    best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]
    
    if(plot){
      # ggplot(aframe, aes(phi, cor)) +
      #   labs(
      #     x = 'phi',
      #     y = 'Correlation'
      #   ) +
      #   geom_vline(xintercept = best_phi, linetype = "dashed") +
      #   geom_point() +
      #   theme_classic()
      aframe$size=1;aframe$size[aframe$phi%in% c(best_cor_phi,best_rmse_phi)]<-2
      ggplot(aframe, aes(rmse, cor,color=phi,size=size)) +
        labs(
          x = 'RSME',
          y = 'Correlation'
        ) +
        geom_point() +
        geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
        theme_classic()
      diffame <-data.frame(diff_corr,diff_rmse,phi=phi_range,size=1)
      diffame$size=1;diffame$size[diffame$phi%in% c(best_cor_phi,best_rmse_phi)]<-2
      ggplot(diffame, aes(diff_rmse, diff_corr,color=phi,size=size)) +
        labs(
          x = 'diff RSME',
          y = 'diff Correlation'
        ) +
        geom_point() +
        geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
        theme_classic()
    }
    
    return(data.frame(best_cor_phi,best_rmse_phi))
  }
  run_reg_lasso_phionly <- function(X, y, scores,
                                    n_folds = 10,
                                    phi_range = seq(0, 1, length = 30)){
    lambda_min <- find_lambda(X, y, plot = F)
    
    # print(paste("Lambda min:", lambda_min))
    print(paste("Lambda 1st SE:", round(lambda_min,4)))
    
    afit <- glmnet(
      X, y, 
      alpha = 1, 
      lambda = lambda_min)
    betas <- afit$beta[,1]
    
    penalties <- scores[match(colnames(X), names(scores))]
    names(penalties) <- colnames(X)
    penalties[is.na(penalties)] <- 0
    penalties <- penalties/1000
    
    asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
    
    correls <- do.call(rbind, lapply(names(asplits), function(x){
      train <- unlist(asplits[setdiff(names(asplits), x)])
      test <- unlist(asplits[x])
      
      do.call(rbind,lapply(phi_range, function(phi){
        lasso_tr <- glmnet(
          X[train,],
          y[train],
          lambda = lambda_min,
          penalty.factor = 1 - penalties * phi)
        
        pred <- predict(lasso_tr,X[test,])
        rmse <- sqrt(apply((y[test]-pred)^2,2,mean))
        if(var(pred) == 0) cor_out <- NA else cor_out <- cor(pred, y[test])
        data.frame(cor= cor_out,rmse=rmse,run=x,phi=phi)
      }))
    }))
    if(any(is.na(correls))) data.frame(best_cor_phi=NA,best_rmse_phi=NA) else find_best_phi(correls, phi_range,plot=F)
  }
  y <- y[, gene]
  scores <- get_scores(gene, ppi)
  out <- run_reg_lasso_phionly(
    X, y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  out$gene<- gene
  return(out)
}
# Define X and y ####
ok_cells <- intersect(
  rownames(rnaseq),
  rownames(kronos))

X <- rnaseq[ok_cells, ]
y <- kronos[ok_cells,]
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))

# save(list=ls(),paste0(wrkfldr,"/wrkspc.RData")
# load(paste0(wrkfldr,"/wrkspc.RData")
# Run command ####


# test phi, corr vs rmse
gene_sample <- c("MYC",sample(colnames(kronos),100))
gene_sample2 <- c("KDM5D","SOX10","FAM50A","RPP25L","PAX8","KRTAP4-11",
                  "EBF1","IRF4","H2BC15","MYB","MDM2","OR4P4"    ,
                  "KRAS","NRAS","HNF1B","OR4C11","EIF1AX","POU2AF1",  
                  "TP63","BRAF","TTC7A","OR4S2")
gene_sample2 <-gene_sample2[gene_sample2%in%colnames(kronos)]

###
library("parallel")

phiout <- mclapply(gene_sample2[1:2],function(gene){
  y <- y[, gene]
  scores <- get_scores(gene, ppi)
  out <- run_reg_lasso_phionly(
    X, y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  out$gene<- gene
  write.csv(out,file=paste0(wrkfldr,gene,'.csv'))
  return(out)
},mc.cores = 10)

do.call(rbind,phiout)
