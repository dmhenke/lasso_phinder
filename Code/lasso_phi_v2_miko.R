# Load R libs ####
library(caret)
library(data.table)
library(ggplot2)
library(glmnet)

## other ####
if(Sys.info()['login']=="Dafydd") setwd("./Code/")


# Source code ####
source("allfunctions.R")


# Load STRING data ####
load("../Data/ppi_w_symbols.RData")
# Loadd Gene infrmation data ####
load("../Data/gene_info.RData")

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
gene_sample6<- c("PRMT5","EGFR","MYC","ERBB2", "CDK4","CDK6","PIK3CA")

res_cnv <- lapply(gene_sample6,function(gene,scorObj='demeter2'){
  save_results_file <- paste0(gene,'_cnv.RData')
  save_results_file_path <-paste0("../Outputs/",scorObj,"/",save_results_file)
  if(save_results_file%in%list.files(paste0("../Outputs/",scorObj,"/"))){
    load(save_results_file_path)
    return(results_cnv)
    }else{
  
  
  if(scorObj=='demeter2'){
    y <- demeter2[, gene]
  } else if(scorObj=='kronos'){
    y <- kronos[, gene]
  }
y <- y[!is.na(y)]

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
tissue <- sample_info[names(y), "lineage"]


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

save(results_cnv,file=save_results_file_path)
return(results_cnv)
}
});names(res_cnv)<-gene_sample6


# plt ####
lapply(names(res_cnv),function(nam,scorObj='demeter2',res_cnv=res_cnv){
 ## prep objs ####
  if(scorObj=='demeter2'){
    y <- demeter2[, nam]
  } else if(scorObj=='kronos'){
    y <- kronos[, nam]
  }
  y <- y[!is.na(y)]
  
  
  
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
  # tissue <- sample_info[names(y), "lineage"]
  # y_resid <- residuals(lm(y ~ tissue))
  # y <- y_resid

  ## Plot ####

  results_cnv <- res_cnv[[nam]]
  tmp <- results_cnv$betas
  tmp$gene <- rownames(tmp)
  tmp$correl <- cor(X_cnv, y)[, 1]
  # reorder
  tmp <- cbind(tmp,gene_info[match(rownames(tmp),gene_info$hgnc_symbol),c("chromosome_name", "start_position")])
  tmp<- tmp[order(tmp$chromosome_name,tmp$start_position),]
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

  # remove X,Y and non cononical chromosomes
  tmp <- tmp[tmp$chromosome_name %in% c(1:22),]
  shaded_chrs2 <- c(as.character(seq(1,22,by=2)),'X',"MT")
  shaded_chrs <- c(as.character(seq(1,22,by=1)),'X',"MT")
  shaded_chrs <- shaded_chrs[shaded_chrs%in% unique(tmp$chromosome_name)]
  shad_block <- do.call(rbind,lapply(shaded_chrs,function(x){
    data.frame(chr=x,xmin=tmp[which(tmp$chromosome_name==x),"rank"][1],xmax=rev(tmp[which(tmp$chromosome_name==x),"rank"])[1])
  }))
  shad_pnt <- rep('gray70',nrow(tmp))
  shad_pnt[which(tmp$chromosome_name  %in% shaded_chrs2)] <- 'black'
  Y_min <- min(tmp$correl)-0.05
  Y_max <- max(tmp$correl)+0.05
  
  out_plt <- ggplot(tmp, aes(rank, correl, label = label)) +
    # annotate("rect", xmin =shad_block$xmin, xmax = shad_block$xmax, ymin = Y_min, ymax = Y_max,alpha = .2)+
    geom_point(alpha=0.5,color=shad_pnt) +
    theme_classic()+
    # ggrepel::geom_label_repel(aes(label=gene,size = abs(betas_pen)),data=glabels,
    #                           box.padding = .3, max.overlaps = Inf, color = "blue2",
    #                           seed = 0,min.segment.length = 0,nudge_y = .1,na.rm = TRUE
    # ) +
    theme(plot.margin=margin(0,0,0,0))+
    # scale_size_continuous("Beta magnitude")+
    scale_size_continuous(NULL)+
    # annotate("text",x=2000,y=-.3,label=paste0("cor=",round(cor.test(aframe[,nam],aframe$y)$estimate,3),
    #                                         "\npval=",round(cor.test(aframe[,nam],aframe$y)$p.value,5)))+
    scale_x_continuous(labels=shaded_chrs, breaks=rowSums(shad_block[,-1])/2)+
    guides(fill="none",color='none',size='none')+
    theme(axis.text.x = element_text(color=rep(c("black","gray70"),length(shaded_chrs)/2)))+
    labs(title=paste(nam,scorObj,"score"),x="Chromosome",y="Correlation")
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"nolab.pdf"),out_plt,width = 12,height = 6)
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"nolab.JPEG"),out_plt,width = 12,height = 6)
  ### Scatter or boxplots explaining primary peak (EGFR amplification is associated with increased dependency)
  # out_plt+ xlim(c(which(tmp$gene=="EGFR")-100,which(tmp$gene=="EGFR")+100))#+ylim()
  # out_plt+  geom_magnify(from = c(1, 0, 5, 1500),to = c(0, 10000, 2.5, 20000),shadow = TRUE)
  # Plot whole chrom
  gene_interest <- nam
  which_interest <- which(tmp$chromosome_name==tmp$chromosome_name[which(tmp$gene==gene_interest)])
  tmpsub <- tmp[which_interest,]
  #only display gene of interest
  glabels_sub<- glabels[which(glabels$gene==nam),]
  out_plt + coord_cartesian(xlim=c(min(which_interest),max(which_interest)),
                            ylim=c(min(tmpsub$correl)-0.05,max(tmpsub$correl)+0.05),
                            expand = F)+
    ggrepel::geom_label_repel(aes(label=gene,size = 1),data=glabels_sub,
                              box.padding = .3, max.overlaps = Inf, color = "blue2",
                              seed = 0,min.segment.length = 0,nudge_y = .05,na.rm = TRUE
    ) +
    labs(title=paste(nam,scorObj,"score\nChr",unique(tmpsub$chromosome_name),"only"))
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"_zoom_",gene_interest,".pdf"),width = 12,height = 6)
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"_zoom_",gene_interest,".JPEG"),width = 12,height = 6)
  ### Zoom into GAB2 region. Show that many genes have similar statistical evidence but lasso identifies GAB2.
  gene_interest <- 'GAB2'
  glabels_sub<- glabels[glabels$gene%in%c(gene_interest,nam),]
  which_interest <- which(tmp$gene==gene_interest)
  tmpsub <- tmp[(which_interest-100):(which_interest+100),]
  out_plt + coord_cartesian(xlim=c(which_interest-100,which_interest+100),
                            ylim=c(min(tmpsub$correl)-0.05,max(tmpsub$correl)+0.05),
                            expand = F)+
    ggrepel::geom_label_repel(aes(label=gene,size = 1),data=glabels_sub,
                              box.padding = .3, max.overlaps = Inf, color = "blue2",
                              seed = 0,min.segment.length = 0,nudge_y = .05,na.rm = TRUE
    ) +
    labs(title=paste(nam,scorObj,"score\n+/- 100 sorrounding genes"))
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"_zoom_",gene_interest,".pdf"),width = 12,height = 6)
  ggsave(paste0("../Outputs/graphics/",nam,"_",scorObj,"_zoom_",gene_interest,".JPEG"),width = 12,height = 6)
  
  
  ## Plot Beta_pen and Beta
  betaP <- ggplot(tmp, aes(rank, betas_pen, label = label)) +
    geom_point(alpha=0.5,color=shad_pnt) +
    theme_classic()+
    ggrepel::geom_label_repel(aes(label=gene),data=glabels,
                              box.padding = .01, max.overlaps = Inf, color = "blue2",
                              seed = 0,min.segment.length = .001,nudge_y = .01,na.rm = TRUE
    ) +
    theme(plot.margin=margin(0,0,0,0))+
    scale_size_continuous(NULL)+
    scale_x_continuous(labels=shaded_chrs, breaks=rowSums(shad_block[,-1])/2)+
    guides(fill="none",color='none',size='none')+
    theme(axis.text.x = element_text(color=rep(c("black","gray70"),length(shaded_chrs)/2)))+
    labs(title=paste(nam,scorObj,"score"),x="Chromosome",y="Correlation")
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
      annotate("text",x=3,y=.2,label=paste0("cor=",round(cor.test(aframe2[,nam],aframe2$y)$estimate,3),
                                                "\npval=",round(cor.test(aframe2[,nam],aframe2$y)$p.value,5)))+
      labs(title="ALK CNV vs CDK6 dependency\npopulation=CDK6 CNV [0.8-1.2]",x=paste("ALK (CNV)"),y=paste("CDK6 (Dependency)"))
    
    
    
    
  }
  # return(out_plt)

})
