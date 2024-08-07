# Load R libs ####
library(biomaRt)
library(caret)
library(data.table)
library(ggplot2)
library(ggpubr)
library(glmnet)


# Source code ####
source("allfunctions.R")


# Load STRING data ####
load("../Data/ppi_w_symbols.RData")


# Load DepMap data ####
load("../Data/global.RData")


# FUNCTIONS ####
plot_manhattan <- function(gene){
  y <- demeter2[, gene]
  y <- y[!is.na(y)]
  
  correl <- cor(
    X_cnv[match(names(y), rownames(X_cnv)), ], y,
    use = "pairwise.complete")[,1]
  
  aframe <- data.frame(
    gene = names(correl),
    gene_info[match(names(correl), gene_info$hgnc_symbol), ],
    correl
  )
  aframe <- aframe[order(aframe$chromosome_name, aframe$start_position),]
  aframe$rank <- 1:nrow(aframe)
  
  subm <- tmp[tmp$omic == "CNV", ]#[tmp$target == gene & tmp$omic == "CNV", ]
  aframe$betas_pen <- subm$betas_pen[match(aframe$gene, subm$gene)]
  
  # aframe <- aframe[!is.na(aframe$chromosome_name), ]
  
  gHattan <- ggplot(aframe,
         aes(rank, correl, color = chromosome_name)) +
    labs(
      title = gene, 
      y = "Pearson correlation coefficient",
      x = "Gene sorted by genomic location"
    ) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    geom_text(
      data = aframe[!is.na(aframe$betas_pen),], 
      aes(label = gene), color = "blue") +
    scale_color_manual(values = rep(c("black", "grey"), 11),guide="none") +
    theme_classic()  
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.png"),plot = gHattan,width = 6,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.pdf"),plot = gHattan,width = 6,height = 4)
  plot(gHattan)
}
plot_manhattan2 <- function(gene="PRMT5",resIn='demeter2_PRMT5_CNV_omic.RData',subplotChr=NA){
  load(paste0("../Outputs/",resIn)) # results_omic
  omic <- strsplit(resIn,"_")[[1]][3]
  y <- demeter2[, gene]
  y <- y[!is.na(y)]

  correl <- results_omic$cor2score
  
  aframe <- data.frame(
    gene = names(correl),
    gene_info[match(names(correl), gene_info$hgnc_symbol), ],
    correl
  )
  aframe <- aframe[order(aframe$chromosome_name, aframe$start_position),]
  aframe$rank <- 1:nrow(aframe)
  
  subm <- results_omic#[resin$omic == "CNV", ]
  subm_betas <- subm$betas[match(aframe$gene, rownames(subm$betas)),]
  subm_betas$betas_pen[which(subm_betas$betas_pen==0)]<- NA;subm_betas$betas[which(subm_betas$betas==0)]<- NA;
  aframe$betas_pen <- subm_betas$betas_pen#[match(aframe$gene, rownames(subm_betas))]
  aframe$betas <- subm_betas$betas#[match(aframe$gene, rownames(subm_betas))]
  aframe$betalogic <- apply(cbind(aframe$betas_pen,aframe$betas),1,function(x){ 
    if(is.na(x[1])&is.na(x[2])) NA else if(!is.na(x[1])&is.na(x[2])) 'beta_pen' else if(!is.na(x[1])&!is.na(x[2])) 'both' else 'beta'
  })
  aframe <- aframe[!is.na(aframe$chromosome_name), ]
  # labs
  lab_x <- "Gene ordered by genomic coordinate"
  lab_y <- paste0("Correlation (r)\n",gene," Dependency (D2) with CNV")
  
  g_full <- ggplot(aframe,
                   aes(rank, correl, color = chromosome_name)) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    # geom_text(
    #   data = aframe[!is.na(aframe$betalogic ),],
    #   aes(label = gene, color = betalogic )) +
    # scale_x_discrete(breaks, labels, limits)+ # consider changing to chr # instead of gene rank
    scale_color_manual(guide='none', breaks=c(1:22,"beta_pen","both","beta"),
                       values = c(rep(c("black", "grey"), 11),"darkblue","purple","darkred")) +
    labs(x=lab_x,y=lab_y)+
    theme_classic()  
  g_fullLab <- g_full +
    geom_text(check_overlap = T,
      data = aframe[which(aframe$betalogic =='beta'),],
      aes(label = gene, color = betalogic , size = abs(betas)))+
    geom_text(check_overlap = T,
      data = aframe[which(!is.na(aframe$betalogic )&aframe$betalogic !='beta'),],
      aes(label = gene, color = betalogic , size = abs(betas_pen)))+
    scale_size(guide = 'none')
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.png"),plot = g_full,width = 10,height = 5)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.pdf"),plot = g_full,width = 10,height = 5)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_labs_demeter2.png"),plot = g_fullLab,width = 10,height = 5)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_labs_demeter2.pdf"),plot = g_fullLab,width = 10,height = 5)
  
  if(is.na(subplotChr)){
  which_chrome <- aframe[which(aframe$gene==gene),"chromosome_name"]
  } else which_chrome <- subplotChr
  which_chromePos <- which(aframe$chromosome_name==which_chrome)
  which_chrome_minX <- min(which_chromePos)
  which_chrome_maxX <- max(which_chromePos)
  
  g_chrGene <- g_fullLab+
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(expand = c(0, 0))+
    labs(x=paste0("Chromosome ",which_chrome))
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_chr",which_chrome,"_demeter2.png"),plot = g_chrGene,width = 2,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_chr",which_chrome,"_demeter2.pdf"),plot = g_chrGene,width = 2,height = 4)
  plot(g_fullLab)
  plot(g_chrGene)
  
  if(gene=="EGFR"&omic=="CNV"){
  # bframe <-  reshape2::melt(aframe,measure.vars=c("betas","betas_pen","correl"))
  bframe <-  aframe;bframe$betas[is.na(bframe$betas)] <- 0;bframe$betas_pen[is.na(bframe$betas_pen)] <- 0;
  plo3b_cor<- ggplot(bframe,aes(x=rank, y=correl)) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    # facet_grid(variable~.)+C
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(expand = c(0, 0))+
    labs(x='',y="Correlation\ncoefficient")+
    theme_classic() 
  # find upper and lower bounds of y <- abs(betas)
  B_hi <- max(c(aframe$betas,aframe$betas_pen),na.rm=T)
  B_low <- min(c(aframe$betas,aframe$betas_pen),na.rm=T)
  B_lim <- max(abs(B_hi),abs(B_low))
  plo3b_b<- ggplot(bframe,aes(x=rank, y=betas)) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(limits = c(-B_lim*1.2,B_lim*1.2),expand = c(0, 0))+
    labs(x='',y="Baseline\ncoefficient")+
    theme_classic() +
    scale_color_manual(guide='none', breaks=c(1:22,"beta_pen","both","beta"),
                       values = c(rep(c("black", "grey"), 11),"darkblue","purple","darkred"))+
    geom_text(check_overlap = T,data = aframe[which(aframe$betalogic =='beta'),],
      aes(label = gene, color = betalogic , size = abs(betas)))+
    geom_text(check_overlap = T,data = aframe[which(!is.na(aframe$betalogic )&aframe$betalogic !='beta'),],
      aes(label = gene, color = betalogic , size = abs(betas_pen)))+theme(legend.position ="none")

  plo3b_bp<- ggplot(bframe,aes(x=rank, y=betas_pen )) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(limits = c(-B_lim*1.2,B_lim*1.2),expand = c(0, 0))+
    labs(x=paste0("Chromosome ",which_chrome),y="Bio-primed\ncoefficient")+
    theme_classic() +
    scale_color_manual(guide='none', breaks=c(1:22,"beta_pen","both","beta"),
                       values = c(rep(c("black", "grey"), 11),"darkblue","purple","darkred"))+
    # geom_text(data = aframe[which(aframe$betalogic =='beta'),],
    #            aes(label = gene, color = betalogic , size = abs(betas)))+
    # geom_text(data = aframe[which(!is.na(aframe$betalogic )&aframe$betalogic !='beta'),],
    #            aes(label = gene, color = betalogic , size = abs(betas_pen)))+theme(legend.position ="none")
    geom_text(check_overlap = T,data = aframe[which(aframe$gene =='GAB2'),],
               aes(label = gene, color = betalogic , size = abs(betas_pen)))+theme(legend.position ="none")
  
  Fig3B <- gridExtra::grid.arrange(plo3b_cor, plo3b_b,plo3b_bp, ncol=1)
  ggsave(filename = paste0("../Outputs/graphics/Fig3B_boxplot_",gene,"_",omic,".pdf"),plot = Fig3B,width = 6,height = 4)
}
}
if(F) {plot_gene <- function(gene,res_path){
  load(res_path)
  subm <- results_omic$betas
  subm <- subm[order(-abs(subm$betas_pen)), ]
  subm <- subm[1:20, ]
  
  subm$gene <- factor(subm$gene, levels = unique(subm$gene))
  
  ggplot(subm, aes(betas_pen, gene, fill = omic)) +
    geom_bar(stat = "identity") +
    labs(
      title = gene, 
      subtitle = "Top 20 multi-omic biomarkers",
      y = "Biomarker",
      x = "Regularized LASSO beta"
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_text(aes(label = gene), size = 3) +
    theme_classic() 
}}
plot_gene <- function(gene = "HNRNPK", target){
  y <- demeter2[, target]
  ok <- intersect(names(y), rownames(cnv))
  aframe <- data.frame(
    cnv = cnv[ok, gene],
    dep = y[ok]
  )
  ggplot(aframe, aes(cnv, dep)) +
    labs(
      x = paste(gene, "copy number"),
      y = paste(target, "dependency [D2]")
    ) +
    geom_smooth(method = lm) +
    geom_point(alpha = 0.75) +
    stat_cor() +
    theme_classic()
}
plot_Xcompare <- function(gene,gene2="GAB2",score_nam="demeter2",X_=cnv,omic=c("CNV","RNA","MUT")[1]){
  if(score_nam=="demeter2") score <- demeter2 else score=kronos
  y <- score[,gene]
  x <- X_[match(names(y),rownames(X_)),]
  
  if(omic=="CNV"){
    omic_cutoff<- 1.2
  } else if(omic=="MUT"){
    omic_cutoff<- 1
  } 

  df <- data.frame(d2=y,log1=x[,gene],log2=x[,gene2]);#colnames(df) <- c("d2",gene,gene2)
  df$d2 <-as.numeric(df$d2)
  df <-df[complete.cases(df),]
  
  breaks <- c('Normal', paste0(gene," only\ngain"), paste0(gene2," only\ngain"), paste0(gene," & ",gene2,"\ngain"))
  breaksVal <- c("FALSE.FALSE","TRUE.FALSE","FALSE.TRUE","TRUE.TRUE")
  names(breaksVal) <- breaks
  breaksVal<-breaksVal[breaksVal%in% unique(interaction(df$log1>=omic_cutoff,df$log2>=omic_cutoff))]
  gbox <-  ggplot(df,aes(y=d2,x=interaction(log1>=omic_cutoff,log2>=omic_cutoff),color=interaction(log1>=omic_cutoff,log2>=omic_cutoff)))+
    geom_jitter(width = .1,size=.5)+
    geom_boxplot(alpha=0.5,outlier.shape = NA)+
    ggpubr::stat_compare_means(comparisons =list(c(1,2),c(1,3),c(2,4)))+
    # scale_x_discrete(labels=c('Normal', paste0(gene," only\ngain"), paste0(gene2," only\ngain"), paste0(gene," & ",gene2,"\ngain")))+
    labs(y=paste0(score_nam,' dependancy: ',gene),x=omic)+
    scale_color_manual(guide='none',values=c("black","blue","darkred","purple3"),
                       breaks=breaksVal)+
                       # labels=c('Normal', paste0(gene," only\ngain"), paste0(gene2," only\ngain"), paste0(gene," & ",gene2,"\ngain")))+
    theme_classic()
  
  ggsave(filename = paste0("../Outputs/graphics/boxplot_",gene,"_",gene2,"_",omic,".png"),plot = gbox,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/boxplot_",gene,"_",gene2,"_",omic,".pdf"),plot = gbox,width = 4,height = 4)
  plot(gbox)
}
plot_betas_multi <- function(res_path="../Outputs/demeter2_PRMT5_results_multiomic.RData"){
  res_info <- strsplit(basename(res_path),'_')[[1]]
  load(res_path)
  plt_betas <-  ggplot(res,aes(x=betas,y=betas_pen,color=omic))+
    geom_hline(yintercept = 0, linetype = 2)+geom_vline(xintercept = 0, linetype = 2)+
    geom_point()+
    ggrepel::geom_text_repel(#point.padding = .5,
                              # min.segment.length = 0,
                              force=1,direction='both',max.overlaps=100,
                              max.time = .3, max.iter = 1e5,#box.padding = 0,
                              data = res[which(abs(res$betas)>=quantile(abs(res$betas),.9) | abs(res$betas_pen)>=quantile(abs(res$betas_pen),.9)),],
                              aes(label = gene, color = omic),size=2) +
    labs(x="LASSO",y="Regularized LASSO",title = paste0(res_info[2]," (",res_info[1],") ~ RNA+CNV+Mutation"))+
    theme_classic()
    ggsave(file= paste0("../Outputs/graphics/BetaScatter_",res_info[2],"_",res_info[1],".png"),plt_betas,height = 8,width = 8)
    ggsave(file= paste0("../Outputs/graphics/BetaScatter_",res_info[2],"_",res_info[1],".pdf"),plt_betas,height = 8,width = 8)
    plot(plt_betas)
}
extractorRData <- function(file, object) {
  #' Function for extracting an object from a .RData file created by R's save() command
  #' Inputs: RData file, object name
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}
# EGFR/GAB2 betas

# Normalize RNA expression data ####
X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(log2(X_cnv))
# offsetlog2 <- max(log2(cnv[which(cnv>2)]))+0.2
# X_cnv[which(cnv<2)] <- cnv[which(cnv<2)] +2^(-offsetlog2)*(2-cnv[which(cnv<2)] )


# Load gene coordinates ####
if(F){mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = colnames(cnv),
  mart = mart)
} else load("~/Projects/lasso_PPI/lasso_ppi/Data/gene_info.RData")
chrs <- as.character(1:22)
gene_info <- gene_info[gene_info$chromosome_name %in% chrs, ]
uniq <- names(which(table(gene_info$hgnc_symbol) == 1))
gene_info <- gene_info[gene_info$hgnc_symbol %in% uniq, ]
gene_info$chromosome_name <- factor(
  gene_info$chromosome_name, levels = chrs)


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut
X_mut <- X_mut[, colSums(X_mut) >= 5]


# Read skewness genes ####
d2_genes <- scan(
  "../Data/skewed_d2_genes.txt",
  what = character(), sep = "\n")

chr_genes <- scan(
  "../Data/skewed_chr_genes.txt",
  what = character(), sep = "\n")


# Load LASSO results ####
load("../Outputs/d2_results_multiomic.RData")
res_d2 <- res
names(res_d2) <- d2_genes

load("../Outputs/chr_results_multiomic.RData")
res_chr <- res
names(res_chr) <- chr_genes


# Assess D2 results ####
res <- res_d2
bad <- which(unlist(lapply(res, class)) == "try-error")
if(length(bad)==0){good <- names(res)}else good <- names(res)[-bad]
d2_combined <- lapply(good, function(x){
  if(is.null(res[[x]])) return(NULL)
  data.frame(target = x, res[[x]])
})
d2_combined <- do.call(rbind, d2_combined)

tmp <- d2_combined
tmp <- tmp[tmp$betas_pen != 0, ]

plot_manhattan(gene = "EGFR")

# Multiomic: PRMT5 ####
# load("../Outputs/demeter2_PRMT5_results_multiomic.RData") # res
load("../Outputs/demeter2_PRMT5_results_multiomic.RData") # res
tmp <- res#[which(res$betas_pen!=0),]
plot_manhattan(gene = "PRMT5")
plot_manhattan2(gene = "PRMT5",resin=res)
plot_betas_multi(res_path="../Outputs/demeter2_PRMT5_results_multiomic.RData")
plot_Xcompare(gene="PRMT5",gene2="MTAP",score_nam="demeter2",X_=cnv,omic="CNV")
plot_Xcompare(gene="PRMT5",gene2="MTAP",score_nam="demeter2",X_=mut,omic="MUT")
# plot_Xcompare(gene="PRMT5",gene2="MTAP",score_nam="demeter2",X=rnaseq,omic='RNA')

# 1) prmt5 boxplot of cnv and d2
#    mtap ...
# 2) scallter copynumber nutral lines (prmt5 & mtap), rna vs d2
# box and scatter plots
lapply(c("PRMT5","MTAP"),function(x){
  samps <- intersect(rownames(cnv),rownames(demeter2))
  dfin <- data.frame(cnv=cnv[samps,x],d2=demeter2[samps,"PRMT5"])
  scat <- ggplot(dfin,aes(x=cnv,y=d2))+
    geom_vline(xintercept = 1,linetype='dashed')+
    annotate("rect",xmin=-.05,xmax=0.05,ymin=-2,ymax=0.3,color='blue',fill='white',linewidth=1.5)+
    geom_point(alpha=0.8)+theme_classic()+
    labs(x=paste0("CNV (",x,")"),y="Demeter2: PRMT5")
  boxUP <- ggplot(dfin,aes(x=cnv>=1.2,y=d2))+geom_boxplot(alpha=0.8)+theme_classic()+
    ggpubr::stat_compare_means(comparisons =list(c(1,2)))+
    labs(#title=paste0("CNV: ",x," vs Demeter2: PRMT5"),
      x=paste0("CNV (",x,") >=1.2"),y="Demeter2: PRMT5")
  boxDOWN <- ggplot(dfin,aes(x=cnv>0.8,y=d2))+geom_boxplot(alpha=0.8)+theme_classic()+
    ggpubr::stat_compare_means(comparisons =list(c(1,2)))+
    labs(      x=paste0("CNV (",x,") >0.8"),y="Demeter2: PRMT5")
  dfin$updwn <- "<=0.8"
  dfin$updwn[ which(dfin$cnv>0.8)]<-  "(0.8-1.2)"
  dfin$updwn[ which(dfin$cnv>=1.2)]<- "1.2"
  dfin$updwn <- factor(dfin$updwn,ordered=T,levels=c("<=0.8","(0.8-1.2)","1.2"))
  boxTRIO <- ggplot(dfin,aes(x=updwn,y=d2))+geom_boxplot(alpha=0.8)+theme_classic()+
    ggpubr::stat_compare_means(comparisons =list(c(1,2),c(2,3)))+
    labs(      x=paste0("CNV (",x,")"),y="Demeter2: PRMT5")
  
  #write graphics
  ggsave(filename = paste0("../Outputs/graphics/BoxP_multiomic_up_",x,"_cnvD2PRMT5.png"),plot = boxUP,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/BoxP_multiomic_down_",x,"_cnvD2PRMT5.png"),plot = boxDOWN,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/BoxP_multiomic_3_",x,"_cnvD2PRMT5.png"),plot = boxTRIO,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/scatter_multiomic_",x,"_cnvD2PRMT5.png"),plot = scat,width = 4,height = 4)
  
  # refine graph of MTAP CNV = 0
  if(x=="MTAP"){
    datsub <- which(dfin$cnv<=0.05)
    df_PRMT5 <- data.frame(cnv=cnv[samps,"PRMT5"],d2=demeter2[samps,"PRMT5"])
    
    SubScat <- ggplot(df_PRMT5[datsub,],aes(x=cnv,y=d2))+
      geom_vline(xintercept = 1,linetype='dashed')+
      ggpubr::stat_cor(label.sep='\n',method='pearson')+
      geom_point(alpha=0.8)+theme_classic()+
      theme(panel.border = element_rect(colour = "blue", fill=NA),axis.line=element_line(color='blue'))+
      labs(x=paste0("CNV (PRMT5)"),y="Demeter2: PRMT5")
    ggsave(filename = paste0("../Outputs/graphics/scatter_CNV_subMTAP_D2PRMT5.png"),plot = SubScat,width = 4,height = 4)
    
  }
  
  # RNA
  df2 <- data.frame(rna=rnaseq[intersect(rownames(rnaseq),rownames(demeter2)),x],d2=demeter2[intersect(rownames(rnaseq),rownames(demeter2)),"PRMT5"])
  # subset to CNV 'neutral"
  df2<- df2[rownames(dfin)[which(dfin$cnv>=0.8 & dfin$cnv<=1.2)],]
  scatR <- ggplot(df2,aes(x=rna,y=d2))+geom_point(alpha=0.8)+theme_classic()+
    labs(title=paste0(x," CNV neutral  lines\nRNA vs Demeter2: PRMT5"),
      x=paste0("RNAseq (",x,")"),y="Demeter2: PRMT5")+
    stat_cor(method = "pearson")#, label.x = 3, label.y = 30)
  ggsave(filename = paste0("../Outputs/graphics/scatter_multiomic_",x,"_rnaD2PRMT5.png"),plot = scatR,width = 4,height = 4)
  
})

plot_manhattan(gene = "MTAP")



plot_manhattan(gene = "MCM2")
plot_manhattan(gene = "RBM39")
plot_manhattan(gene = "DDX46")

# e.g. EGFR 
load("~/Projects/lasso_PPI/lasso_ppi/Outputs/demeter2_EGFR_results_CNV_omic.RData")
# d2_combined <- res
tmp <- d2_combined <- res
# tmp <- tmp[tmp$betas_pen != 0, ]

# EGFR CNV only
load("~/Projects/lasso_PPI/lasso_ppi/Outputs/demeter2/EGFR_cnv.RData") # results_cnv
bets <- re$betas
ok_cell <- intersect(rownames(demeter2),rownames(X_cnv))
cor_D2cnv <-  cor(
  X_cnv[ok_cell, ], demeter2[ok_cell,'EGFR'],
  use = "pairwise.complete")[,1]
bets$cor <- cor_D2cnv[match(rownames(bets),names(cor_D2cnv))]
bets$gene <- rownames(bets)
bets <-cbind(bets,gene_info[match(bets$gene,gene_info$hgnc_symbol),])

chrqry <- bets$chromosome_name[which(bets$hgnc_symbol==gene2)]
ggplot(bets[which(bets$chromosome_name==chrqry),],aes(x=start_position, y=cor))+
  geom_point()
bets11 <-bets[which(bets$chromosome_name==chrqry),]
bets11m <- reshape2::melt(bets11,measure.vars=c("betas","betas_pen","cor"))
bets11m$variable <- factor(bets11m$variable,levels=c("cor","betas","betas_pen"),labels = c("Correlation","LASSO","Regularized LASOO"))
Gcor_betas <- ggplot(bets11m,aes(x=start_position, y=value))+
  geom_point(size=.5)+theme_classic()+
  labs(x=paste0('Position on Chr ',chrqry),y=NULL)+
  facet_wrap(~variable,ncol=1)

ggsave(filename = paste0("../Outputs/graphics/StackCorBetas_",gene,"_",gene2,".png"),plot = Gcor_betas,width = 4,height = 4)
ggsave(filename = paste0("../Outputs/graphics/StackCorBetas_",gene,"_",gene2,".pdf"),plot = Gcor_betas,width = 4,height = 4)








# EGFR ####
{
# plot_gene(gene = "EGFR")
plot_manhattan2(gene = "EGFR",resIn = "../Outputs/demeter2_EGFR_CNV_omic.RData")

d2_EGFR<- demeter2[,'EGFR']
egfr_cnv <- X_cnv[,"EGFR"]
egfr_cells <- intersect(names(d2_EGFR),names(egfr_cnv))
d2_EGFR<-d2_EGFR[egfr_cells]
egfr_cnv<-egfr_cnv[egfr_cells]
ggplot(data.frame(d2=d2_EGFR,x=egfr_cnv),aes(x=x,y=d2))+geom_point()+
  geom_smooth(method = "lm",se=F) +
  geom_point() +
  ggpubr::stat_cor(label.sep='\n') +
  theme_classic() +
  labs(x="EGRF [log2(CNV)]",y="EGFR dependency [D2]")
ggsave(filename = paste0("../Outputs/FigS4_CNV_EGFRvsD2.pdf"),height=4,width=4)

}
# DMAP1 (figures 4) ####
{
## Fig 4A ####
plot_manhattan2(gene = "DMAP1",resIn = "../Outputs/demeter2_DMAP1_CNV_omic.RData")
## Fig 4B ####
  for(gene in 'DMAP1'){
  res <- extractorRData("../Outputs/demeter2_DMAP1_CNV_omic.RData",'results_omic')
  basB <- rownames(res$betas)[which(res$betas$betas!=0)]
  bioB <- rownames(res$betas)[which(res$betas$betas_pen!=0)]
  basB <-basB[!basB%in%gene]
  bioB <-bioB[!bioB%in%gene]
  # boxplot
  # all hits
  bioB <- bioB[bioB%in%colnames(demeter2)]
  basB <- basB[basB%in%colnames(demeter2)]

  cor_pen <- cor(
    demeter2[, bioB], demeter2[, gene],
    use = "pairwise.complete")[,1]
  cor_reg <- cor(
    demeter2[, basB], demeter2[, gene],
    use = "pairwise.complete")[,1]
  aframe <- data.frame(
    class = c(rep("Bio-Primed", length(cor_pen)),
              rep("Baseline", length(cor_reg))),
    cor = c(cor_pen, cor_reg),
    gene=c(names(cor_pen),names(cor_reg)))
  ggplot(aframe, aes(x=class,cor, color = class)) +
    labs(
      x = "Co-dependency with target",
      y = "Percentile",
      color="Process"
      ) +
    geom_hline(yintercept = 0,linetype='dashed')+
    geom_boxplot() +
    theme_classic()
  ggsave(filename = paste0("../Outputs/graphics/FigS4C_boxplt_codep_DMAP1vsbiomarkers.pdf"),height=4,width=4)
  # boxplot
  ggplot(aframe, aes(class, cor, label = gene)) +
    labs(
      y = "Correlation (r)\nDependency DMAP1 [D2] vs biomarker [D2]",
      x = "Biomarkers derived from process",
      color="Process") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.1) +
    ggrepel::geom_text_repel(max.overlaps = 100, aes(color = class)) +
    scale_color_manual(values = c("red", "blue")) +
    theme_classic()
  ggsave(filename = paste0("../Outputs/graphics/FigS4B_codep_DMAP1vsbiomarkers.pdf"),height=4,width=4)
  
  # unique to lasso
  bioB <- bioB[!bioB%in%basB][bioB[!bioB%in%basB]%in%colnames(demeter2)]
  basB <- basB[!basB%in%bioB][basB[!basB%in%bioB]%in%colnames(demeter2)]
  
  cor_pen <- cor(
    demeter2[, bioB], demeter2[, gene],
    use = "pairwise.complete")[,1]
  cor_reg <- cor(
    demeter2[, basB], demeter2[, gene],
    use = "pairwise.complete")[,1]
  
  # # wilcox rank sum test
  # wilcox.test(cor_pen, cor_reg)
  
  aframe <- data.frame(
    class = c(rep("Bio-Primed", length(cor_pen)),
              rep("Baseline", length(cor_reg))),
    cor = c(cor_pen, cor_reg),
    gene=c(names(cor_pen),names(cor_reg)))
  # ECDF
  ggplot(aframe, aes(cor, color = class)) +
    labs(
      x = "Co-dependency with target",
      y = "Percentile",
      color="Unique gene\nto process"
    ) +
    stat_ecdf() +
    xlim(-0.3, 0.3) +
    theme_classic()
  ggsave(filename = paste0("../Outputs/graphics/Fig4D_boxplt_codep_DMAP1vsbiomarkers.pdf"),height=4,width=4)
  # boxplot
  ggplot(aframe, aes(class, cor, label = gene)) +
    labs(
      y = "Correlation (r)\nDependency DMAP1 [D2] vs biomarker [D2]",
      x = "Biomarkers derived from process",
      color="Process"
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.1) +
    ggrepel::geom_text_repel(max.overlaps = 100, aes(color = class)) +
    scale_color_manual(values = c("red", "blue")) +
    theme_classic()
  ggsave(filename = paste0("../Outputs/graphics/Fig4B_codep_DMAP1vsbiomarkers.pdf"),height=4,width=4)
  
  }
## Fig C ####
 # against BRD8  CNV and D2 
 # against MCRS1 CNV and D2
lapply(c("BRD8","MCRS1"),function(x){
  d2_DMAP1<- demeter2[,'DMAP1']
  # CNV
  x_cnv <- X_cnv[,x]
  x_cells <- intersect(names(d2_DMAP1),names(x_cnv))
  d2_DMAP1 <- d2_DMAP1[x_cells]
  x_cnv <-    x_cnv[x_cells]
  ggplot(data.frame(d2=d2_DMAP1,x=x_cnv),aes(x=x,y=d2))+geom_point()+
    geom_smooth(method = "lm",se=F) +
    geom_point(alpha=0.75) +
    ggpubr::stat_cor(label.sep='\n') +
    theme_classic() +
    labs(x=paste0(x," [log2(CNV)]"),y="DMAP1 dependency [D2]")
  ggsave(filename = paste0("../Outputs/graphics/Fig4C_scatt_DMAP1_D2vs",x,"_CNV.pdf"),height=4,width=4)
  # D2
  d2_x <- demeter2[,x]
  d2_DMAP1<- demeter2[,'DMAP1']
  ggplot(data.frame(d2=d2_DMAP1,x=d2_x),aes(x=x,y=d2))+geom_point()+
    geom_smooth(method = "lm",se=F) +
    geom_point(alpha=0.75) +
    ggpubr::stat_cor(label.sep='\n') +
    theme_classic() +
    labs(x=paste0(x," dependency [D2]"),y="")
  ggsave(filename = paste0("../Outputs/graphics/Fig4C_scatt_DMAP1_D2vs",x,"_D2.pdf"),height=4,width=4)
})

## Fig 4D ####


}



plot_gene("AURKA")
plot_gene("TUBB")


gene <- "AURKA"
subm <- tmp[tmp$target == gene, ]
n_mut <- rowSums(
  X_mut[, intersect(subm$gene, colnames(X_mut))])

ok <- intersect(rownames(X_mut), rownames(demeter2))
boxplot(split(demeter2[ok, gene], n_mut[ok]))


# Assess Chronos results ####
res <- res_chr
bad <- which(unlist(lapply(res, class)) == "try-error")
good <- names(res)[-bad]
chr_combined <- lapply(good, function(x){
  if(is.null(res[[x]])) return(NULL)
  data.frame(target = x, res[[x]])
})
chr_combined <- do.call(rbind, chr_combined)

tmp <- chr_combined
tmp <- tmp[tmp$betas_pen != 0, ]

