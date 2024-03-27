# All functions for project #
if(!'wrkfldr'%in%ls()) wrkfldr < getwd()

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

## Given matrix of correlation between predicted and observed values from
# cross-validation across range of phi values, find best phi ####
find_best_phi <- function(correls, phi_range, plot = T){
  median_correls <- unlist(lapply(split(correls$cor,correls$phi), median))
  median_rmse <- unlist(lapply(split(correls$rmse,correls$phi), median))
  
  aframe <- data.frame(
    phi = phi_range, 
    cor = median_correls,
    rmse = median_rmse)
  
  afit <- lm(cor ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- predict(afit, aframe)
  diff_corr <- aframe$cor - preds
  afit2 <- lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds2 <- predict(afit2, aframe)
  diff_rmse <- aframe$rmse - preds2
  
  best_cor_phi <- phi_range[which(diff_corr == max(diff_corr))]
  best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]
  
  if(plot){
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

## Run lasso ####
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
  # sniff test #
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
  if(length(unique(dim(correls) == dim(na.omit(correls)))) == 2){
    print("Missing values in correlation.")
    return(NA)
  }
  
  best_phi <- find_best_phi(correls, phi_range)
  best_phi_cor <- best_phi$best_cor_phi
  print(paste("Phi best:", round(best_phi_cor,4)))
  
  # Update lambda and rerun LASSO
  afit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi_cor)
  betas_pen <- afit$beta[,1]
  
  # write outs ####
  # phis ##
  write.csv(best_phi,file=paste0(wrkfldr,gene,'_phibest.csv'))
  # cor & rmse & phi ##
  write.csv(correls,file=paste0(wrkfldr,gene,'_phibest.csv'))
  # betas ##
  write.csv(data.frame(betas, betas_pen),file=paste0(wrkfldr,gene,'_phibest.csv'))
  
  return(list(best_phi,
              best_phi_cor,
              correlations = correls,
              betas = data.frame(betas, betas_pen)))
}


