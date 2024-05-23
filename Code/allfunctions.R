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
