# R load R libs
library(glmnet)

# Given a single gene symbol extract PPI confidence scores
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

# Given X and y find best lambda using standard glmnet implementation
find_lambda <- function(X, y, plot = T){
  fitcv <- cv.glmnet(
    X, y,
    alpha = 1,
    lambda = NULL)
  if(plot) plot(fitcv, xvar = "lambda", label = T)
  fitcv$lambda.min
}

# Given matrix of correlation between predicted and observed values from
# cross-validation across range of phi values, find best phi
find_best_phi <- function(correls, phi_range, plot = T){
  median_correls <- apply(correls, 2, median)

  aframe <- data.frame(
    phi = phi_range,
    cor = median_correls)

  afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- predict(afit, aframe)
  diffs <- aframe$cor - preds
  best_cor <- median_correls[which(diffs == max(diffs))]
  best_phi <- phi_range[which(diffs == max(diffs))]

  if(plot){
    ggplot(aframe, aes(phi, cor)) +
      labs(
        x = 'phi',
        y = 'Correlation'
      ) +
      geom_vline(xintercept = best_phi, linetype = "dashed") +
      geom_point() +
      theme_classic()
  }

  return(c(best_phi, best_cor))
}

# Given X, y and scores, run both baseline and final lasso models
run_reg_lasso <- function(X, y, scores,
                          n_folds = 10,
                          phi_range = seq(0, 1, length = 30)){
  
  lambda_min <- find_lambda(X, y, plot = F)
  
  print(paste("Lambda min:", lambda_min))
  
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
  
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  
  correls <- do.call(rbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    
    unlist(lapply(phi_range, function(phi){
      lasso_tr <- glmnet(
        X[train,],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
        pred <- predict(
          lasso_tr,
          X[test,])
      
        cor(pred, y[test])
      }))
    }))
  
    if(length(unique(dim(correls) == dim(na.omit(correls)))) == 2){
      print("Missing values in correlation.")
      return(NA)
    }
  
    tmp <- find_best_phi(correls, phi_range, plot = F)
    best_phi <- tmp[1]
    best_phi_correl <- tmp[2]
  
    print(
      paste("Phi best:", best_phi,
            "with correlation:", best_phi_correl))
  
    afit <- glmnet(
          X,
          y,
          alpha = 1,
          lambda = lambda_min,
          penalty.factor = 1 - penalties * best_phi)
    betas_pen <- afit$beta[,1]
  
    return(
      list(best_phi,
           best_phi_correl,
           correlations = correls,
           betas = data.frame(betas, betas_pen)))
  }
