# Load R libs ####
library(biomaRt)
library(caret)
library(data.table)
library(ggplot2)
library(ggpubr)
library(glmnet)


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


# Load LASSO results ####
files <- list.files("../Outputs/", full.names = T)
files <- files[grep("demeter2", files)]

res <- lapply(files, function(x){
  load(x)
  data.frame(
    target = strsplit(basename(x), "_", fixed = T)[[1]][2],
    res
  )
})
names(res) <- unlist(lapply(files, function(x)
  strsplit(basename(x), "_", fixed = T)[[1]][2]))

cnv_markers <- lapply(res, function(x)
  x[x$omic == "CNV" & x$betas_pen != 0, ])
cnv_markers <- cnv_markers[unlist(lapply(cnv_markers, nrow)) > 0]

freq_markers <- unlist(lapply(cnv_markers, nrow))
max_cor <- unlist(lapply(cnv_markers, function(x)
  max(abs(x$correl))))

freq_markers <- data.frame(
  gene = names(freq_markers),
  freq = freq_markers,
  max_cor)

ggplot(freq_markers, aes(x = freq)) +
  labs(x = "Number of CNV biomarkers") +
  geom_histogram() +
  theme_classic()


# Check FLT4 ####
plot_top_hits <- function(gene){
  tmp <- res[[gene]]
  tmp <- tmp[tmp$betas_pen != 0, ]
  tmp$marker <- paste(tmp$gene, tmp$omic, sep = "_")
  tmp$marker <- factor(
    tmp$marker,
    levels = tmp$marker[order(tmp$betas_pen)])
  tmp <- tmp[order(tmp$marker), ]
  tmp <- tmp[c(1:20, (nrow(tmp) - 10): nrow(tmp)), ]
  
  ggplot(tmp, aes(betas_pen, marker, fill = omic)) +
    labs(title = gene) +
    geom_bar(stat = "identity") +
    labs(x = "Beta coefficient", y = "Biomarker") +
    theme_classic()  
}
plot_top_hits("HPRT1")

plot_manhattan <- function(gene){
  y <- demeter2[, gene]
  y <- y[!is.na(y)]
  ok <- intersect(names(y), rownames(X_cnv))
  correl <- cor(X_cnv[ok, ], y[ok])[,1]
  correl <- data.frame(
    gene = names(correl),
    correl = correl)
  
  tmp <- res[[gene]]
  tmp <- tmp[tmp$betas_pen != 0, ]
  tmp <- tmp[tmp$omic == "CNV", ]
  
  correl$betas_pen <- tmp$betas_pen[match(correl$gene, tmp$gene)]
  
  correl <- correl[match(gene_info$hgnc_symbol, correl$gene), ]
  correl$pos <- gene_info$start_position
  correl$chr <- gene_info$chromosome_name
  
  correl <- correl[order(correl$chr, correl$pos), ]
  
  correl$rank <- 1:nrow(correl)
  
  ggplot(correl, aes(rank, correl)) +
    labs(
      title = gene,
      x = "Genes sorted by genomic location",
      y = "Pearson correlation"
    ) +
    geom_point() +
    geom_text(data = correl[correl$betas_pen != 0, ],
              aes(label = gene, size = abs(betas_pen)),
              color = "red") +
    theme_classic() +
    theme(legend.position = "none")
}
plot_manhattan("FLT4")

