# Load R libs ####
library(caret)
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

# Define outcome ####
gene <- "MYC"
y <- kronos[, gene]
y <- y[!is.na(y)]
# Define working set of samples ###
ok_cells <- intersect(
  rownames(X),
  names(y)
)

X <- X[ok_cells, ]
y <- y[ok_cells]
# Generate scores ###
scores <- get_scores(gene, ppi)

# Run LASSO ####
### LOAD RUN MYC RNA RESULTS ####
file_myc_rnaKRONOS <- "kronos_MYC_RNA_omic.RData"
if(file_myc_rnaKRONOS%in% list.files("../Outputs/")){
  load(paste0("../Outputs/",file_myc_rnaKRONOS))
} else{
results <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))
save(results,file=paste0("../Outputs/",file_myc_rnaKRONOS))
}

# Show how phi was inferred: deviation from RMSE ####
# 1) rmse over phi test range at 10x cross validation
# 2) linear fit of median RMSE
# 3) select phi with greatest deviation from median RMSE
tmp <- results$correlations; rownames(tmp) <- seq(0, 1, length = 30)
# asplit <- split(1:nrow(tmp), tmp$run)
# tmp <- do.call(cbind, lapply(asplit, function(x) tmp$cor[x]))
# tmp <- apply(tmp, 2, function(x)
#   (x - mean(x))/sd(x))

phi_range <- seq(0, 1, length = 30)
median_cor <- apply(tmp, 1, median)

aframe <- data.frame(
  phi = phi_range, 
  cor = median_cor)

fitline<- coef(lm(cor ~ phi, data = aframe[c(1,30),]))

ggplot(aframe, aes(phi, cor)) +
  labs(
    title = "MYC [Chronos] vs RNA expression",
    y = "Median standardized RMSE", #standardized between x-validation by selected phi
    x = "phi"
  ) +
  geom_abline(intercept=fitline[1],slope=fitline[2],color='red') +
  geom_point() +
  geom_vline(xintercept = results[[1]],linetype='dashed') +
  theme_classic()
ggsave(filename = paste0("../Outputs/Fig1X_myc_phiselect.png"),height=4,width=4)
ggsave(filename = paste0("../Outputs/Fig1X_myc_phiselect.pdf"),height=4,width=4)

# Show non-zero coefficients ####
aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  cor = cor(X, y)[,1])

aframe$gene <- factor(
  aframe$gene, 
  levels = aframe$gene[order(aframe$betas_pen)])
ggplot(aframe[order(abs(aframe$betas_pen),decreasing = T)[1:40], ], aes(betas_pen, gene)) +
  labs(
    y = "Informative & relevant features",
    x = "Regularized LASSO coefficient"
  ) +
  geom_bar(stat = "identity") +
  theme_classic()
ggsave(filename = paste0("../Outputs/Fig2A_myc_top40betaPen_barplz.pdf"),height=4,width=4)


# Compare to regular correlation coefficients ####
aframe$label <- aframe$gene
aframe$label[aframe$betas_pen == 0] <- NA

# myc beta results
myc_xmin <- -0.05
myc_xmax <- 0.05
myc_ymin <- -0.05
myc_ymax <- 0.05

g_betaPen_myc <- ggplot(aframe, aes(betas, betas_pen,
                   label = label)) +
  labs(
    y = "Regularized LASSO beta coefficient",
    x = "LASSO beta coefficient"
  ) +
  annotate('rect',xmin=myc_xmin,xmax=myc_xmax,ymin=myc_ymin,ymax=myc_ymax,color='black',alpha=.9,fill='white',linetype="dashed")+
  geom_point() +
  # ggrepel::geom_label_repel(max.overlaps = 20) +
  ggrepel::geom_label_repel(data = aframe[which(aframe$gene=="MYC"),],
                            aes(label = gene))+
  theme_classic()
ggsave(filename = paste0("../Outputs/Fig2B_myc_betaVSbetapen.pdf"),height=4,width=4,plot = g_betaPen_myc)

#zoom in on non-myc results
g_betaPen_mycSub<- ggplot(aframe[aframe$gene != "MYC", ], aes(betas, betas_pen, 
                                           color = cor,
                                           label = label)) +
  geom_vline(xintercept = 0,linetype='dashed')+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(
    y = "Regularized LASSO coefficient",
    x = "LASSO coefficient"
  ) +
  scale_color_gradient2(
    low = "blue", mid = "grey", high = "red",
    breaks = seq(-0.5, 0.5, length = 11)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  theme_classic()
ggsave(filename = paste0("../Outputs/Fig2C_myc_betaVScorrSub.pdf"),height=6,width=6,plot = g_betaPen_mycSub)



plot_gene <- function(gene, y){
  subm <- data.frame(
    gene = X[, gene],
    y
  )  
  
  ggplot(subm, aes(gene, y)) +
    labs(
      x = paste(gene, "RNA levels"),
      y = "MYC dependency [Chronos]"
    ) +
    geom_smooth(method = "lm",se=F) +
    geom_point() +
    ggpubr::stat_cor() +
    theme_classic() 
}
g_rnaMYC<-plot_gene("MYC", y)
g_rnUBR5<-plot_gene("UBR5", y)
g_rnaWEE1<-plot_gene("WEE1", y)
ggsave(filename = paste0("../Outputs/Fig2D_RNA_MYC.pdf"),height=4,width=4,plot = g_rnaMYC)
ggsave(filename = paste0("../Outputs/Fig2D_RNA_UBR5.pdf"),height=4,width=4,plot = g_rnUBR5)
ggsave(filename = paste0("../Outputs/Fig2D_RNA_WEE1.pdf"),height=4,width=4,plot = g_rnaWEE1)

# Run one more time ####
results_new <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

aframe <- data.frame(
  results$betas,
  results_new$betas
)
g_betPen_1v2<- ggplot(aframe,
       aes(betas_pen, betas_pen.1)) +
  labs(
    title = "Regularized betas",
    x = "Run 1",
    y = "Run 2"
  ) +
  geom_abline(intercept = 0,slope=1)+
  geom_point(alpha=0.5) +
  ggpubr::stat_cor() +
  theme_classic()
ggsave(filename = paste0("../Outputs/FigS2A_betaPen_run1vs2.pdf"),height=4,width=4,plot = g_betPen_1v2)


# Run LASSO with MYC score set to 0 ####
scores["MYC"] <- 0

set.seed(202407)
results_myc0 <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  results_myc0$betas
)

aframe$label <- aframe$gene
aframe$label[aframe$betas_pen == 0] <- NA

g_regBeta_set0 <- ggplot(aframe,
       aes(betas_pen, betas_pen.1, label = label)) +
  labs(
    title = "Regularized betas",
    x = "Original",
    y = "MYC score set to 0"
  ) +
  geom_abline(intercept = 0,slope=1)+
  geom_point(alpha=0.5) +
  ggrepel::geom_label_repel() +
  ggpubr::stat_cor() +
  theme_classic()
ggsave(filename = paste0("../Outputs/FigS3_betaPen_runvsZero.pdf"),height=4,width=4,plot = g_regBeta_set0)

# Compare biomarkers to predict MYC dependency in D2 ####
myc_d2 <- demeter2[, "MYC"]
myc_d2 <- myc_d2[!is.na(myc_d2)]

X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))

ok_cells <- intersect(
  rownames(X),
  names(myc_d2)
)

X <- X[ok_cells, ]
myc_d2 <- myc_d2[ok_cells]


# D2results <- run_reg_lasso(
#   X, myc_d2, scores,
#   n_folds = 10, phi_range = seq(0, 1, length = 30))
# d2cor <- apply(D2results$correlations, 1, median)
cor_ChRNA <-  cor(
  X, y[ok_cells],
  use = "pairwise.complete")[,1]
res_chroMYC <-aframe
res_chroMYC$cor <- cor_ChRNA


genes_bio <- res_chroMYC$gene[res_chroMYC$betas_pen != 0]
genes_bet <- res_chroMYC$gene[res_chroMYC$betas != 0]
genes_cor <- tail(res_chroMYC$gene[order(abs(res_chroMYC$cor))], max(length(genes_bio),length(genes_bet)))

library(ggVennDiagram)
ggVennDiagram(list(biolasso=genes_bio, lasso=genes_bet,corr=genes_cor))
intersect(genes_bio, genes_cor)

fitControl <- trainControl(
  method = "cv", number = 10, savePredictions = T)
# unique biolasso
set_bio <- setdiff(setdiff(genes_bio, genes_bet),genes_cor)#setdiff(genes_bio, genes_cor)#
subm <- data.frame(myc_d2, X[, set_bio])
ngene_biolasso <- length(set_bio)
afit_bio <- train(
  myc_d2 ~ ., data = subm, method = "lm", trControl = fitControl)
# unique lasso
set_lasso <- setdiff(setdiff(genes_bet,genes_bio),genes_cor)#setdiff(genes_bet, genes_cor)#
subm <- data.frame(myc_d2, X[, set_lasso])
ngene_lasso <- length(set_lasso)
afit_las <- train(
  myc_d2 ~ ., data = subm, method = "lm", trControl = fitControl)
# unique corr
set_cor <-setdiff(setdiff(genes_cor,genes_bio),genes_bet)#setdiff(genes_cor,genes_bio)#
subm <- data.frame(myc_d2, X[, set_cor])
ngene_corr <- length(set_cor)
afit_cor <- train(
  myc_d2 ~ ., data = subm, method = "lm", trControl = fitControl)




aframe <- data.frame(
  biolas = afit_bio$resample$Rsquared,
  lasso = afit_las$resample$Rsquared,
  cor = afit_cor$resample$Rsquared,
  rep=1:10);colnames(aframe) <- c(paste0('biolas\n',ngene_biolasso),
                                  paste0('lasso\n',ngene_lasso),paste0('corr\n',ngene_corr),'rep')
aframe <- reshape2::melt(aframe,'rep')

ggplot(aframe, aes(x=variable, y=value)) +
  labs(
    title = "Predicting MYC dependency in D2",
    y = "R-squared",
    x = "Biomarkers from"
  ) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.1) +
  ggpubr::stat_compare_means(method = "t.test",comparisons =list(c(1,3),c(2,3))) +
  theme_classic()

