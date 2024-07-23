# All functions for project #

# R libs
library(ggplot2)
library(glmnet)

# Define helper functions ####
## Get Association Scores ####
get_scores <- function(gene, ppi){
  tmp <- ppi[(ppi$gene1 %in% c(gene) | 
                ppi$gene2 %in% c(gene)), ]
  tmp <- tmp[tmp$gene1 != "" & tmp$gene2 != "", ]
  tmp <- na.omit(tmp)
  scores <- as.numeric(tmp[["combined_score"]][tmp$gene1 == gene])
  names(scores) <- as.character(tmp[["gene2"]][tmp$gene1 == gene])
  scores[gene] <- 1000
  return(scores)
}

## Calculate best lambda ####
find_lambda <- function(X, y, plot = F){
  fitcv <- cv.glmnet(
    X, y, 
    alpha = 1, 
    lambda = NULL)
  if(plot) plot(fitcv, xvar = "lambda", label = T)
  
  fitcv$lambda.min # more conservative vs 'lambda.min'
}

## Scale input data (by gene)
scale <- function(matr){
  apply(matr, 2, function(x)
    (x - mean(x))/sd(x))
}

## Given matrix of correlation between predicted and observed values from

## cross-validation across range of phi values, find best phi ####
find_best_phi_rmse <- function(correls, phi_range, plot = F){
  median_rmse <- apply(correls, 1, median)
  
  aframe <- data.frame(
    phi = phi_range, 
    rmse = median_rmse)
  
  afit <- lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- predict(afit, aframe)
  diff_rmse <- aframe$rmse - preds
  
  best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]

  
  if(plot){
    ggplot(aframe, aes(phi, rmse)) +
      labs(
        x = 'Phi',
        y = 'RMSE'
      ) +
      geom_point() +
      geom_vline(xintercept = best_rmse_phi) +
      theme_classic()
    # rmse vs phi
    ggplot(aframe, aes(phi, cor,color=phi,size=size)) +
      labs(
        x = 'Phi',
        y = 'Correlation'
      ) +
      geom_point() +
      geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
      theme_classic()+geom_vline(xintercept = best_cor_phi)
    ggplot(aframe, aes(phi, rmse,color=phi,size=size)) +
      labs(
        x = 'Phi',
        y = 'RMSE'
      ) +
      geom_point() +
      geom_abline(intercept =afit$coefficients[1],slope = afit$coefficients[2] )+
      theme_classic()+geom_vline(xintercept = best_cor_phi)
  }
  
  return(best_rmse_phi)
}

## Get rmse values from cross-validation ####
get_rmse_from_xvalid <- function(
    X, y, penalties, lambda_min, phi_range, n_folds = 10){
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  
  correls <- do.call(cbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    
    do.call(rbind, lapply(phi_range, function(phi){
      lasso_tr <- glmnet(
        X[train,],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
      
      pred <- predict(lasso_tr, X[test,])
      rmse <- sqrt(apply((y[test] - pred)^2, 2, mean))
      
      rmse
    }))
  }))
  correls <- apply(correls, 2, function(x)
    (x - mean(x))/sd(x))
  
  return(correls)
}

## Run lasso ####
run_reg_lasso <- function(X, y, scores,
                          n_folds = 10,
                          phi_range = seq(0, 1, length = 30)){
  lambda_min <- find_lambda(X, y, plot = F)
  
  print(paste("Lambda min:", round(lambda_min,4)))
  
  afit <- glmnet(
    X, y, 
    alpha = 1, 
    lambda = lambda_min)
  betas <- afit$beta[,1]
  
  if(sum(betas) == 0){
    print("All betas are zero.")
    return(NA)
  }
  
  penalties <- scores[match(colnames(X), names(scores))]
  names(penalties) <- colnames(X)
  penalties[is.na(penalties)] <- 0
  penalties <- penalties/max(scores)
  
  correls <- get_rmse_from_xvalid(
    X, y, penalties, phi_range = phi_range, lambda_min = lambda_min, n_folds = n_folds)
  
  if(length(unique(dim(correls) == dim(na.omit(correls)))) == 2){
    print("Missing values in correlation.")
    return(NA)
  }
  

  best_phi <- find_best_phi_rmse(correls, phi_range, plot = F)

  print(paste("Best phi based on RMSE:", round(best_phi, 4)))
  
  # Update lambda and rerun LASSO
  afit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi)
  betas_pen <- afit$beta[,1]
  
  return(list(best_phi,
              correlations = correls,
              betas = data.frame(betas, betas_pen)))
}

# Run multiomic for list of genes on D2 ####
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
  
  # All 'omics in one lasso
  colnames(X_rna_ok) <- paste0(colnames(X_rna_ok),"_rna")
  colnames(X_cnv_ok) <- paste0(colnames(X_cnv_ok),"_cnv")
  colnames(X_mut_ok) <- paste0(colnames(X_mut_ok),"_mut")
  cat_omics <- cbind(scale(X_rna_ok),scale(X_cnv_ok),scale(X_mut_ok))
  system.time({
  results_omics <- run_reg_lasso(
    cat_omics, y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  })
  
  
  
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
  
  genes <- unlist(lapply(betas, function(x) colnames(x)))
  omic <- unlist(lapply(names(betas), function(x){rep(x, ncol(betas[[x]])) }))
  
  
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

# Run singleomic for list of genes on D2 ####
run_analysisSingle <- function(whichY="demeter2", gene,omics=list(X_rna,X_cnv,X_mut),which_omic=c("RNA","CNV","Mut",NULL)[4]){
  if(whichY=="demeter2")y <- demeter2[,gene] else y <- kronos[, gene]
  
  y <-na.omit(y)
  lapply(omics,function(omic){
    min_omic <- min(omic[,sample(colnames(omic),5)])
    max_omic <- max(omic[,sample(colnames(omic),5)])
    # Find overlapping cell lines ####
  ok_cells <- intersect(names(y), rownames(omic))

  # Remove features without variance ####
    omic_OK  <- omic[ok_cells, ]
    if(is.null(which_omic)){
      which_omic <-'tmp'}
    if(min_omic==0 & max_omic==1|tolower(which_omic)=="mut"){
      # "Mutation"
      omic_OK <- omic_OK[, colSums(omic_OK) >= 5]
      which_omic<-"Mut"
    } else if(min_omic>=0 & max_omic>1|tolower(which_omic)=="rna"){
      # "RNA"
      omic_OK <- omic_OK[, apply(omic_OK, 2, var) > 0]
      which_omic <-"RNA"
    }else  if (min_omic<0 & max_omic>0|tolower(which_omic)=="cnv"){
      # "CNV"
      omic_OK <- omic_OK[, apply(omic_OK, 2, var) > 0]
      which_omic<- "CNV"
    }
    
    y_ok <- y[ok_cells]
    # check if analysis already run
    output_file <- paste0("../Outputs/",whichY,"_",gene,"_",which_omic,"_omic.RData")
    if(!basename(output_file) %in% list.files("../Outputs/")){
      
  # Run LASSO ####
  scores <- get_scores(gene, ppi)
  
  results_omic <- run_reg_lasso(
    scale(omic_OK), y_ok, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  results_omic$cor2score <- cor(
    omic_OK, y_ok,
    use = "pairwise.complete")[,1]
  results_omic$omic <- which_omic
  # write out resutls
  save(results_omic,file=output_file)
  return(results_omic)
    }
})
}
