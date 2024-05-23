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


# Normalize RNA expression data ####
X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X_rna <- X


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(log2(X_cnv))


# Load gene coordinates ####
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = colnames(cnv),
  mart = mart)
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
  
  subm <- tmp[tmp$target == gene & tmp$omic == "CNV", ]
  aframe$betas_pen <- subm$betas_pen[match(aframe$gene, subm$gene)]
  
  aframe <- aframe[!is.na(aframe$chromosome_name), ]
  
  ggplot(aframe,
         aes(rank, correl, color = chromosome_name)) +
    labs(
      title = gene, 
      y = "Pearson correlation coefficient",
      x = "Gene sorted by genomic location"
    ) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point() +
    geom_label(
      data = aframe[!is.na(aframe$betas_pen),], 
      aes(label = gene), color = "blue") +
    scale_color_manual(values = rep(c("black", "grey"), 11)) +
    theme_classic()  
}
plot_manhattan(gene = "MCM2")
plot_manhattan(gene = "RBM39")
plot_manhattan(gene = "DDX46")


plot_gene <- function(gene){
  subm <- tmp[tmp$target == gene, ]
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

