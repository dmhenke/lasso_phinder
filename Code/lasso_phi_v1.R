# Load R libs ####
library(biomaRt)
library(data.table)
library(ggplot2)
library(glmnet)


# Load String data ####
# ppi <- fread("../Data/9606.protein.links.v12.0.txt")
# 
# ppi$protein1 <- gsub("9606.", "", ppi$protein1, fixed = T)
# ppi$protein2 <- gsub("9606.", "", ppi$protein2, fixed = T)
# 
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_info <- getBM(
#   attributes = c("ensembl_peptide_id", "hgnc_symbol","chromosome_name","start"), mart = mart)
# nam_chr <- c(as.character(c(1:22)),"X","Y","MT")
# gene_info$chromosome_name <- factor(gene_info$chromosome_name , ordered = TRUE, 
#        levels = c(nam_chr,sort(unique(gene_info$chromosome_name)[!unique(gene_info$chromosome_name)%in%nam_chr])))
# save(gene_info,file="../Data/gene_info.RData")
# 
# ppi$gene1 <- gene_info$hgnc_symbol[match(ppi$protein1, gene_info$ensembl_peptide_id)]
# ppi$gene2 <- gene_info$hgnc_symbol[match(ppi$protein2, gene_info$ensembl_peptide_id)]
# save(ppi,file="../Data/ppi_w_symbols.RData")

load("../Data/ppi_w_symbols.RData")
load("../Data/gene_info.RData")

# Load data: expressions & scores ####
# load("C:/Users/Dafydd/Documents/Projects/lasso_PPI/global.RData")
load(paste0("../Data/global.RData"))


# Load functions ####
source("./allfunctions.R")


# Define X and y ####
# Dependency score source
# Y
scor_src <- 'demeter2'
if(scor_src == 'demeter2'){
  dep_scor <- demeter2} else if(scor_src == 'kronos'){
    dep_scor <- kronos}
# X
exp_src <- 'cnv'
if(exp_src=='rnaseq'){
  ind_exp <- rnaseq
}else if(exp_src=='cnv'){
  ind_exp <- na.omit(cnv)
}

ind_exp <- na.omit(ind_exp)


ok_cells <- intersect(
  rownames(ind_exp),
  rownames(dep_scor))

X_ind_exp <- ind_exp[ok_cells, ]
y_dep_scor <- dep_scor[ok_cells,]
if(exp_src=='rnaseq'){
  expressed <- apply(X_ind_exp, 2, function(x) mean(x > 0))
  X_ind_exp <- X_ind_exp[, expressed > 0.95]
  X_ind_exp <- apply(X_ind_exp, 2, function(x)
    (x - mean(x))/sd(x))
}

# test phi, corr vs rmse
gene_sample <- c("MYC",sample(colnames(dep_scor),100))
gene_sample2 <- c("KDM5D","SOX10","FAM50A","RPP25L","PAX8","KRTAP4-11",
                  "EBF1","IRF4","H2BC15","MYB","MDM2","OR4P4"    ,
                  "KRAS","NRAS","HNF1B","OR4C11","EIF1AX","POU2AF1",  
                  "TP63","BRAF","TTC7A","OR4S2")
gene_sample <-gene_sample[gene_sample%in%colnames(dep_scor)]
gene_sample2 <-gene_sample2[gene_sample2%in%colnames(dep_scor)]

# BRCA1 assocated genes
gene_sample3 <-c("PALB2","BARD1","BRIP1","RAD50","RAD51C","RAD51D","RAD54L","ATM","ATR","ATRX","FANCONI","FANCA", "FANCB", "FANCC","CHEK1", "CHEK2", "BLM","NBN","MRE11A")
gene_sample3 <-gene_sample3[gene_sample3%in%colnames(dep_scor)]

gene_sample4 <- c("MYC","SF3B1","BRCA1","BRCA2")
# gene_sample5<- c("KRAS","TP53","RASGRF2", "NRAS",""
gene_sample6<- c("EGFR","MYC","ERBB2", "CDK4/6","PIK3CA")
                                  


library("parallel")

phiout <- mclapply(gene_sample4,function(gene){
  save_results_file <- paste0(gene,'.RData')
  if(save_results_file%in%list.files(paste0("../Outputs/",scor_src,"/"))){
    print(paste('alredy processed',gene))
    NULL
  }else{
    print(paste('Processing',gene))
    y <- na.omit(y_dep_scor[, gene])
    X <- X_ind_exp[match(names(y),rownames(X_ind_exp)), ]
    scores <- get_scores(gene, ppi)
    out <- run_reg_lasso(
      X, y, scores,
      n_folds = 10, phi_range = seq(0, 1, length = 30))
    save(out,file=paste0("../Outputs/",scor_src,"/",gene,'.RData'))
    return(out)
  }
},mc.cores = 10)




# PLOTTING ####


tmp <- results_cnv$betas
tmp$gene <- rownames(tmp)
tmp$correl <- cor(X_cnv, y)[, 1]
tmp$rank <- 1:nrow(tmp)

tmp$label <- tmp$gene
tmp$label[tmp$betas_pen == 0] <- NA
tmp$label[tmp$gene == "SF3B1"] <- "SF3B1"


glabels <- tmp[order(abs(tmp$betas_pen),decreasing = T)[1:20],]


ggplot(tmp, aes(rank, correl, label = label)) +
  geom_point() +
  # geom_label(aes(size = abs(betas_pen))) +
  theme_classic()+
  ggrepel::geom_label_repel(aes(label=gene,size = abs(betas_pen)),data=glabels,
                            box.padding = .3, max.overlaps = Inf, color = "blue2",
                            seed = 0,min.segment.length = 0,nudge_y = .1,na.rm = TRUE
  ) +
  labs(title="SF3B1",x="Gene order",y="Correlation")




# BELOW HERE IS TEMP ####
# Load R libs ####
library(caret)
library(data.table)
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
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(X_cnv)


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut


# Define outcome ####
gene_sample6<- c("EGFR","MYC","ERBB2", "CDK4","CDK6","PIK3CA")

res_cnv <- lapply(gene_sample6[5],function(gene){
  y <- demeter2[, gene]
  # y <- kronos[, gene]
  y <- y[!is.na(y)]
  
  tissue <- sample_info[names(y), "lineage"]
  #y_resid <- residuals(lm(y ~ tissue))
  #y <- y_resid
  
  
  # Find overlapping cell lines ####
  ok_cells <- intersect(names(y), rownames(X_rna))
  ok_cells <- intersect(ok_cells, rownames(X_cnv))
  ok_cells <- intersect(ok_cells, rownames(X_mut))
  
  
  # Remove features without variance ####
  X_rna  <- X_rna[ok_cells, ]
  X_rna <- X_rna[, apply(X_rna, 2, var) > 0]
  
  X_mut <- X_mut[ok_cells, ]
  X_mut <- X_mut[, colSums(X_mut) >= 5]
  
  X_cnv  <- X_cnv[ok_cells, ]
  X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]
  
  y <- y[ok_cells]
  
  
  # Run LASSO ####
  scores <- get_scores(gene, ppi)
  
  # results_rna <- run_reg_lasso(
  #   X_rna, y, scores,
  #   n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  results_cnv <- run_reg_lasso(
    X_cnv, y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  # results_mut <- run_reg_lasso(
  #   X_mut, y, scores,
  #   n_folds = 10, phi_range = seq(0, 1, length = 30))
  
  
  return(results_cnv)
  
});names(res_cnv)<-gene_sample6

# plt ####
lapply(names(res_cnv),function(nam){
  ## prep objs ####
  y <- demeter2[, nam]
  # y <- kronos[, gene]
  y <- y[!is.na(y)]
  
  tissue <- sample_info[names(y), "lineage"]
  #y_resid <- residuals(lm(y ~ tissue))
  #y <- y_resid
  
  
  # Find overlapping cell lines ####
  ok_cells <- intersect(names(y), rownames(X_rna))
  ok_cells <- intersect(ok_cells, rownames(X_cnv))
  ok_cells <- intersect(ok_cells, rownames(X_mut))
  
  
  # Remove features without variance ####
  X_rna  <- X_rna[ok_cells, ]
  X_rna <- X_rna[, apply(X_rna, 2, var) > 0]
  
  X_mut <- X_mut[ok_cells, ]
  X_mut <- X_mut[, colSums(X_mut) >= 5]
  
  X_cnv  <- X_cnv[ok_cells, ]
  X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]
  
  y <- y[ok_cells]
  
  ## Plot ####
  
  results_cnv <- res_cnv[[nam]]
  tmp <- results_cnv$betas
  tmp$gene <- rownames(tmp)
  tmp$correl <- cor(X_cnv, y)[, 1]
  tmp$rank <- 1:nrow(tmp)
  
  tmp$label <- tmp$gene
  tmp$label[tmp$betas_pen == 0] <- NA
  tmp$label[tmp$gene == nam] <- nam
  
  
  glabels <- tmp[order(abs(tmp$betas_pen),decreasing = T)[1:20],]
  # corr
  aframe <- data.frame(
    X_cnv[, c(nam, "KRAS", "NRAS","ALK")],
    y
  )
  
  out_plt <- ggplot(tmp, aes(rank, correl, label = label)) +
    geom_point() +
    # geom_label(aes(size = abs(betas_pen))) +
    theme_classic()+
    ggrepel::geom_label_repel(aes(label=gene,size = abs(betas_pen)),data=glabels,
                              box.padding = .3, max.overlaps = Inf, color = "blue2",
                              seed = 0,min.segment.length = 0,nudge_y = .1,na.rm = TRUE
    ) +
    annotate("text",x=2000,y=-.3,label=paste0("cor=",round(cor.test(aframe[,nam],aframe$y)$estimate,3),
                                              "\npval=",round(cor.test(aframe[,nam],aframe$y)$p.value,5)))+
    labs(title=nam,x="Gene order",y="Correlation")
  if(scatterp){
    aframe$CDK6_norm <- F
    aframe$CDK6_norm[which(aframe$CDK6>=0.8&aframe$CDK6<=1.2)]<-T
    ggplot(aframe,aes(y=y,x=aframe[,"ALK"],color=CDK6_norm))+
      theme_classic()+      geom_hline(yintercept =0)+
      geom_vline(xintercept =1)+geom_point()+
      geom_smooth(method="lm",se=F)+
      # scale_x_continuous(breaks = c(.25,.5,1,2,3,4),labels =c("1/4",'1/2','1','2','3','4'),trans="log2")+
      labs(title="ALK CNV vs CDK6 dependency",x=paste("ALK (CNV)"),y=paste("CDK6 (Dependency)"))
    
    # plot only normal 
    aframe2 <- aframe[which(aframe$CDK6_norm),]
    ggplot(aframe2,aes(y=y,x=aframe2[,"ALK"],color=CDK6_norm))+
      theme_classic()+      
      geom_hline(yintercept =0)+
      geom_vline(xintercept =1)+geom_point()+
      geom_smooth(method="lm",se=F)+
      annotate("text",x=3,y=.2,label=paste0("cor=",round(cor.test(aframe2[,"ALK"],aframe2$y)$estimate,3),
                                            "\npval=",round(cor.test(aframe2[,"ALK"],aframe2$y)$p.value,5)))+
      labs(title="ALK CNV vs CDK6 dependency\npopulation=CDK6 CNV [0.8-1.2]",x=paste("ALK (CNV)"),y=paste("CDK6 (Dependency)"))
    
    
    ggplot(aframe, aes(interaction(ALK > 1.2, CDK6 > 1.2), y)) + geom_boxplot()
    aframe$ALK_norm<-F; aframe$ALK_norm[which(aframe$ALK>=0.8&aframe$ALK<=1.2)]<-T
    ggplot(aframe,aes(x=y,color=ALK_norm,fill=CDK6_norm,group=paste(ALK_norm,CDK6_norm)))+geom_histogram(position = 'identity',alpha=0.5)+
      geom_vline(aes(xintercept=mean(y)))+
      theme_classic()+
      facet_grid(paste(CDK6_norm,ALK_norm)~.)
    
  }
  # return(out_plt)
  ggsave(paste0("../Outputs/graphics/",nam,"_kronos.pdf"),out_plt,width = 8,height = 4)
})
