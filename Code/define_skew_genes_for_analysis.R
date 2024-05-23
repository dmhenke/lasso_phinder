# Load R libs ####
library(moments)


# Load Depmap data ####
load("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/global.RData")


# Define good genes for D2 ####
tmp <- demeter2
nas <- apply(tmp, 2, function(x) mean(is.na(x)))
tmp <- tmp[, nas < 0.1]

num_dependent <- apply(tmp, 2, function(x) mean(x < -0.5, na.rm = TRUE))
tmp <- tmp[, num_dependent > 0.05]

skew <- apply(tmp, 2, skewness, na.rm = TRUE)
good <- names(which(abs(skew) < 0.5))

write(
  good,
  file = "/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/skewed_d2_genes.txt", 
  sep = "\n")


# Define good genes for Chronos ####
tmp <- kronos

num_dependent <- apply(tmp, 2, function(x) mean(x < -0.5, na.rm = TRUE))
tmp <- tmp[, num_dependent > 0.05]

skew <- apply(tmp, 2, skewness, na.rm = TRUE)
good <- names(which(abs(skew) < 0.5))

write(
  good,
  file = "/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/skewed_chr_genes.txt", 
  sep = "\n")
