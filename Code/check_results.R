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
    geom_label(
      data = aframe[!is.na(aframe$betas_pen),], 
      aes(label = gene), color = "blue") +
    scale_color_manual(values = rep(c("black", "grey"), 11),guide="none") +
    theme_classic()  
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.png"),plot = gHattan,width = 6,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.pdf"),plot = gHattan,width = 6,height = 4)
  
}
plot_manhattan2 <- function(gene,resin=res){
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
  
  subm <- resin[resin$omic == "CNV", ]
  subm$betas_pen[which(subm$betas_pen==0)]<- NA;subm$betas[which(subm$betas==0)]<- NA;
  aframe$betas_pen <- subm$betas_pen[match(aframe$gene, subm$gene)]
  aframe$betas <- subm$betas[match(aframe$gene, subm$gene)]
  aframe$bestlogic <- apply(cbind(aframe$betas_pen,aframe$betas),1,function(x){ 
    if(is.na(x[1])&is.na(x[2])) NA else if(!is.na(x[1])&is.na(x[2])) 'beta_pen' else if(!is.na(x[1])&!is.na(x[2])) 'both' else 'beta'
  })
  aframe <- aframe[!is.na(aframe$chromosome_name), ]
  # aframe$chromosome_color <- "grey";aframe$chromosome_color[aframe$chromosome_name%in%seq(1,22,by=2)] <- "black";
  # labs
  lab_x <- "Gene ordered by genomic coordinate"
  lab_y <- paste0("Correlation (r)\n",gene," Dependency (D2) with CNV")
  
  g_full <- ggplot(aframe,
                   aes(rank, correl, color = chromosome_name)) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    # geom_label(
    #   data = aframe[!is.na(aframe$bestlogic),],
    #   aes(label = gene, color = bestlogic)) +
    scale_color_manual(guide='none', breaks=c(1:22,"beta_pen","both","beta"),
                       values = c(rep(c("black", "grey"), 11),"darkblue","purple","darkred")) +
    labs(x=lab_x,y=lab_y)+
    theme_classic()  
  g_fullLab <- g_full +
    geom_label(
      data = aframe[!is.na(aframe$bestlogic),],
      aes(label = gene, color = bestlogic)) 
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.png"),plot = g_full,width = 6,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_demeter2.pdf"),plot = g_full,width = 6,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_labs_demeter2.png"),plot = g_fullLab,width = 6,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_labs_demeter2.pdf"),plot = g_fullLab,width = 6,height = 4)
  
  

  which_chrome <- aframe[which(aframe$gene==gene),"chromosome_name"]
  which_chromePos <- which(aframe$chromosome_name==which_chrome)
  which_chrome_minX <- min(which_chromePos)
  which_chrome_maxX <- max(which_chromePos)
  
  g_chrGene <- g_full+
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(expand = c(0, 0))+
    labs(x=paste0("Chromosome ",which_chrome))
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_chr",which_chrome,"_demeter2.png"),plot = g_chrGene,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/Manhattan_",gene,"_chr",which_chrome,"_demeter2.pdf"),plot = g_chrGene,width = 4,height = 4)
  
}
plot_gene <- function(gene){
  subm <- tmp#[tmp$target == gene, ]
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
}
plot_Xcompare <- function(gene,gene2="GAB2",score_nam="demeter2",X_cnv=cnv){
  if(score_nam=="demeter2") score <- demeter2 else score=kronos
  y <- score[,gene]
  x <- X_cnv[match(names(y),rownames(X_cnv)),]
  # norm_gene <- x[,gene]>0.8&x[,gene]<1.2
  # norm_gene1 <- x[,gene2]>0.8&x[,gene2]<1.2
  # 
  # norm_gene <- x[,gene]<=0.8
  # norm_gene1 <- x[,gene2]<=0.8
  
  # norm_gene <- x[,gene]>=1.2
  # norm_gene1 <- x[,gene2]>=1.2
  # 
  # df <- as.data.frame(rbind(cbind(d2=y,log1=norm_gene,Gene=gene),cbind(d2=y,log1=norm_gene1,Gene=gene2)))
  # df$d2 <-as.numeric(df$d2)
  # df <-df[!is.na(df$log1),]
  # ggplot(df,aes(y=d2,x=paste(Gene,log1)))+geom_boxplot()+
  #   facet_wrap(~Gene)+
  #   ggpubr::stat_compare_means()
  
  df <- data.frame(d2=y,log1=x[,gene],log2=x[,gene2]);colnames(df) <- c("d2",gene,gene2)
  df$d2 <-as.numeric(df$d2)
  df <-df[complete.cases(df),]
  gbox <- ggplot(df,aes(y=d2,x=interaction(EGFR>=1.2,GAB2>=1.2)))+
    geom_jitter(width = .1,size=.5)+
    geom_boxplot(alpha=0.5,outlier.shape = NA)+
    ggpubr::stat_compare_means(comparisons =list(c(1,2),c(1,3),c(2,4)))+
    scale_x_discrete(labels=c('Normal', paste0(gene," only\ngain"), paste0(gene2," only\ngain"), paste0(gene," & ",gene2,"\ngain")))+
    labs(y=paste0(score_nam,' dependancy: ',gene),x="CNV")+
    theme_classic()
  
  ggsave(filename = paste0("../Outputs/graphics/boxplot_",gene,"_",gene2,".png"),plot = gbox,width = 4,height = 4)
  ggsave(filename = paste0("../Outputs/graphics/boxplot_",gene,"_",gene2,".pdf"),plot = gbox,width = 4,height = 4)
  
  }
# EGFR/GAB2 betas

# Normalize RNA expression data ####
X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
offsetlog2 <- max(log2(cnv[which(cnv>2)]))+0.2
X_cnv[which(cnv<2)] <- cnv[which(cnv<2)] +2^(-offsetlog2)*(2-cnv[which(cnv<2)] )
X_cnv <- na.omit(log2(X_cnv))


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
good <- names(res)[-bad]
d2_combined <- lapply(good, function(x){
  if(is.null(res[[x]])) return(NULL)
  data.frame(target = x, res[[x]])
})
d2_combined <- do.call(rbind, d2_combined)

tmp <- d2_combined
tmp <- tmp[tmp$betas_pen != 0, ]



plot_manhattan(gene = "EGFR")

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









plot_gene(gene = "EGFR")


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

