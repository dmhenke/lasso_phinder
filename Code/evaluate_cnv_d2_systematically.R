# Load R libs ####
library(biomaRt)
library(data.table)
library(ggplot2)
library(glmnet)
library(parallel)


# Set working directory on taco ####
if(Sys.info()[[6]]=="Dafydd"){ setwd('../Code/')}else setwd("/storage/thinc/git_repos/lasso_phinder/Code")


# Source code ####
source("allfunctions.R")


# Load DepMap data ####
load("../Data/global.RData")


# Get genome coordinates ####
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = colnames(cnv),
  mart = mart)

chrs <- as.character(1:22)
gene_info <- gene_info[gene_info$chromosome_name %in% chrs, ]

uniq <- names(which(table(gene_info$hgnc_symbol) == 1))
gene_info <- gene_info[which(gene_info$hgnc_symbol %in% uniq), ]
gene_info$chromosome_name <- factor(gene_info$chromosome_name, levels = chrs)


# Define mutation table ####
tmp <- fread("../Data/Damaging_Mutations.csv")
mut <- data.matrix(tmp[,-1])
mut <- apply(mut, 2, function(x)
  as.numeric(x > 0))
rownames(mut) <- tmp[[1]]
X_mut <- mut
X_mut <- X_mut[, colSums(X_mut) >= 5]


# Define genes ####
genes <- list.files("../Outputs/cnv_d2/")
genes <- gsub(".RData", "", fixed = T, genes)


# Define functions ####
get_effect <- function(gene){
  load(paste0(
    "/storage/thinc/git_repos/lasso_phinder/Outputs/cnv_d2/", gene,".RData"))
  
  y <- demeter2[, gene]
  y <- y[!is.na(y)]
  
  tmp <- res$betas
  tmp$gene <- rownames(tmp)
  
  genes_pen <- tmp[tmp$betas_pen != 0, ]
  genes_pen <- genes_pen[order(-abs(genes_pen$betas_pen)), ]
  genes_pen <- genes_pen$gene[genes_pen$betas_pen > 0]
  genes_pen <- setdiff(genes_pen, gene)
  if(length(genes_pen) > 10) genes_pen <- genes_pen[1:10] 
  
  genes_reg <- tmp[tmp$betas != 0, ]
  genes_reg <- genes_reg[order(-abs(genes_reg$betas)), ]
  genes_reg <- genes_reg$gene[genes_reg$betas > 0]
  genes_reg <- setdiff(genes_reg, gene)
  if(length(genes_reg) > 10) genes_reg <- genes_reg[1:10]
  
  common <- intersect(genes_pen, genes_reg)
  genes_pen <- setdiff(genes_pen, common)
  genes_reg <- setdiff(genes_reg, common)
  
  ok <- intersect(names(y), rownames(mut))
  y <- y[ok]
  
  mut_small <- mut[ok, ]
  mut_small <- mut_small[, colSums(mut_small) > 0]
  
  genes_pen <- intersect(genes_pen, colnames(mut_small))
  genes_reg <- intersect(genes_reg, colnames(mut_small))
  
  treat_pen <- rowSums(mut_small[, genes_pen]) > 0
  treat_reg <- rowSums(mut_small[, genes_reg]) > 0
  
  coef_pen <- c(NA, NA)
  coef_reg <- c(NA, NA)
  
  if(sum(treat_pen) >= 5){
    coef_pen <- coefficients(summary(lm(y ~ treat_pen)))[2, c(1, 4)]
  }
  if(sum(treat_reg) >= 5){
    coef_reg <- coefficients(summary(lm(y ~ treat_reg)))[2, c(1, 4)]
  } 
  
  aframe <- rbind(
    data.frame(
      class = "penalized",
      num_mut = sum(treat_pen),
      coef = coef_pen[1],
      pval = coef_pen[2]),
    data.frame(
      class = "regular",
      num_mut = sum(treat_reg),
      coef = coef_reg[1],
      pval = coef_reg[2])
  )
  
  return(aframe)
}
plot_manhattan <- function(gene){
  load(paste0(
    "/storage/thinc/git_repos/lasso_phinder/Outputs/cnv_d2/", gene,".RData"))
  
  betas <- data.frame(
    gene = rownames(res$betas),
    res$betas)
  
  y <- demeter2[, gene]
  y <- y[!is.na(y)]
  
  ok <- intersect(names(y), rownames(cnv))
  
  tmp <- res$betas
  tmp$gene <- rownames(res$betas)
  
  tmp$correl <- cor(cnv[ok, tmp$gene], y[ok])[,1]
  tmp$chr <- gene_info$chromosome_name[match(tmp$gene, gene_info$hgnc_symbol)]
  tmp$pos <- gene_info$start_position[match(tmp$gene, gene_info$hgnc_symbol)]
  tmp <- tmp[order(tmp$chr, tmp$pos), ]
  tmp$rank <- 1:nrow(tmp)
  
  ggplot(tmp[!is.na(tmp$chr), ],
         aes(rank, correl, color = chr)) +
    labs(
      title = gene,
      y = "Correlation coefficient",
      x = "Genes sorted by genomic location"
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    scale_color_manual(values = rep(c("black", "grey"), 11)) +
    ggrepel::geom_text_repel(
      max.overlaps = 100,
      data = tmp[tmp$betas_pen != 0 & !is.na(tmp$chr), ],
      aes(label = gene, size = abs(betas_pen)), color = "red") +
    ggrepel::geom_text_repel(
      max.overlaps = 100,
      data = tmp[tmp$betas != 0 & !is.na(tmp$chr), ],
      aes(label = gene, size = abs(betas)), color = "blue") +
    theme_classic() +
    theme(legend.position = "none")
}
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
get_correlation <- function(gene){
  load(paste0(
    "/storage/thinc/git_repos/lasso_phinder/Outputs/cnv_d2/", gene,".RData"))
  
  tmp <- res$betas
  tmp$gene <- rownames(tmp)
  
  genes_pen <- tmp[tmp$betas_pen != 0, ]
  genes_pen <- genes_pen[order(-abs(genes_pen$betas_pen)), ]
  genes_pen <- genes_pen$gene[genes_pen$betas_pen > 0]
  genes_pen <- setdiff(genes_pen, gene)
  if(length(genes_pen) > 20) genes_pen <- genes_pen[1:20]
  
  genes_reg <- tmp[tmp$betas != 0, ]
  genes_reg <- genes_reg[order(-abs(genes_reg$betas)), ]
  genes_reg <- genes_reg$gene[genes_reg$betas > 0]
  genes_reg <- setdiff(genes_reg, gene)
  if(length(genes_reg) > 20) genes_reg <- genes_reg[1:20]
  
  common <- intersect(genes_pen, genes_reg)
  genes_pen <- setdiff(genes_pen, common)
  genes_reg <- setdiff(genes_reg, common)
  
  genes_pen <- intersect(genes_pen, colnames(demeter2))
  genes_reg <- intersect(genes_reg, colnames(demeter2))
  
  cor_pen <- cor(
    demeter2[, genes_pen], demeter2[, gene],
    use = "pairwise.complete")[,1]
  cor_reg <- cor(
    demeter2[, genes_reg], demeter2[, gene],
    use = "pairwise.complete")[,1]
  
  list(cor_pen, cor_reg)  
}
plot_codep <- function(gene1 = "DMAP1", gene2 = "CDK11B"){
  aframe <- data.frame(
    gene1 = demeter2[, gene1],
    gene2 = demeter2[, gene2])
  
  ggplot(aframe, aes(gene1, gene2)) +
    labs(
      x = paste(gene1, "dependency [D2]"),
      y = paste(gene2, "dependency [D2]")
    ) +
    geom_smooth(method = lm) +
    geom_point(alpha = 0.75) +
    stat_cor() +
    theme_classic()
}


# Run mutation analysis ####
res <- lapply(genes, function(x)
  try(get_effect(x)))
names(res) <- genes
res <- res[which(unlist(lapply(res, class)) != "try-error")]

merged <- do.call(rbind, res)
merged$gene <- unlist(lapply(rownames(merged), function(x)
  strsplit(x, ".", fixed = T)[[1]][1]))

tmp_pen <- merged[merged$class == "penalized", ]
tmp_pen$gene <- names(res)

tmp_reg <- merged[merged$class == "regular", ]
tmp_reg$gene <- names(res)

tmp_combined <- data.frame(tmp_pen, tmp_reg)


# Create mutation volcano plots ####
ggplot(merged, aes(coef, -log10(pval), 
                   color = pval < 0.01)) +
  facet_wrap(~ class) +
  scale_color_manual(values = c("grey", "red")) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = merged[merged$pval < 0.01 & !is.na(merged$pval), ],
    aes(label = gene)
  ) +
  theme_classic()

table(merged$class, merged$pval < 0.01)

ggplot(merged, aes(pval, color = class)) +
  geom_density() +
  theme_classic()

tmp_combined <- data.frame(tmp_pen, tmp_reg)
ggplot(tmp_combined, aes(-log10(pval), -log10(pval.1))) +
  geom_point() +
  theme_classic()


# Highlight mutation examples ####
gene <- "DMAP1"
load(paste0(
  "/storage/thinc/git_repos/lasso_phinder/Outputs/cnv_d2/", gene,".RData"))

betas <- data.frame(
  gene = rownames(res$betas),
  res$betas)

y <- demeter2[, gene]
y <- y[!is.na(y)]

ok <- intersect(names(y), rownames(cnv))

tmp <- res$betas
tmp$gene <- rownames(res$betas)

tmp$correl <- cor(cnv[ok, tmp$gene], y[ok])[,1]
tmp$chr <- gene_info$chromosome_name[match(tmp$gene, gene_info$hgnc_symbol)]
tmp$pos <- gene_info$start_position[match(tmp$gene, gene_info$hgnc_symbol)]
tmp <- tmp[order(tmp$chr, tmp$pos), ]
tmp$rank <- 1:nrow(tmp)

genes_pen <- tmp[tmp$betas_pen != 0, ]
genes_pen <- genes_pen[order(-abs(genes_pen$betas_pen)), ]
genes_pen <- genes_pen$gene[genes_pen$betas_pen > 0]
genes_pen <- setdiff(genes_pen, gene)
if(length(genes_pen) > 10) genes_pen <- genes_pen[1:10] 

genes_reg <- tmp[tmp$betas != 0, ]
genes_reg <- genes_reg[order(-abs(genes_reg$betas)), ]
genes_reg <- genes_reg$gene[genes_reg$betas > 0]
genes_reg <- setdiff(genes_reg, gene)
if(length(genes_reg) > 10) genes_reg <- genes_reg[1:10]

common <- intersect(genes_pen, genes_reg)
genes_pen <- setdiff(genes_pen, common)
genes_reg <- setdiff(genes_reg, common)

ok <- intersect(names(y), rownames(mut))
y <- y[ok]

mut_small <- mut[ok, ]
mut_small <- mut_small[, colSums(mut_small) > 0]

genes_pen <- intersect(genes_pen, colnames(mut_small))
genes_reg <- intersect(genes_reg, colnames(mut_small))

treat_pen <- rowSums(mut_small[, genes_pen]) > 0
treat_reg <- rowSums(mut_small[, genes_reg]) > 0

subm <- data.frame(
  y, treat_pen, treat_reg
)
subm <- reshape2::melt(subm, id.var = "y")
ggplot(subm, aes(value, y)) +
  labs(
    y = paste(gene, "dependency [D2]"),
    x = "Mutation in >0 biomarkers"
  ) +
  facet_wrap(~ variable) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.1) +
  stat_compare_means() +
  theme_classic()


# Run co-dependency analysis ####
correls <- lapply(genes, function(x)
  try(get_correlation(x)))
names(correls) <- genes
correls <- correls[which(unlist(lapply(correls, class)) != "try-error")]

cor_pen <- unlist(lapply(correls, function(x) x[[1]]))
cor_reg <- unlist(lapply(correls, function(x) x[[2]]))

aframe <- data.frame(
  class = c(rep("penalized", length(cor_pen)),
            rep("regular", length(cor_reg))),
  cor = c(cor_pen, cor_reg))

wilcox.test(cor_pen, cor_reg)

ggplot(aframe, aes(cor, color = class)) +
  labs(
    x = "Co-dependency with target",
    y = "Percentile"
  ) +
  stat_ecdf() +
  xlim(-0.3, 0.3) +
  theme_classic()

ggsave("../Outputs/codependency_ecdf.pdf", height = 4, width = 4)

pvals <- lapply(correls, function(x) wilcox.test(x[[1]], x[[2]])$p.value)
head(sort(unlist(pvals)))

plot_manhattan("DMAP1")

ggsave("../Outputs/dmap1_manhattan.pdf", height = 6, width = 12)

subm <- correls[["DMAP1"]]
subm <- data.frame(class = c(
  rep("penalized", length(subm[[1]])),
  rep("standard", length(subm[[2]]))),
  cor = unlist(subm),
  gene = c(names(subm[[1]]), names(subm[[2]])))

ggplot(subm, aes(class, cor, label = gene)) +
  labs(
    y = "Correlation with DMAP1 dependency",
    x = "Biomarkers derived from"
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.1) +
  ggrepel::geom_text_repel(max.overlaps = 100, aes(color = class)) +
  scale_color_manual(values = c("red", "blue")) +
  theme_classic()

ggsave("../Outputs/dmap1_boxplot.pdf", height = 6, width = 6)


p1 <- plot_gene("MCRS1", "DMAP1")
p2 <- plot_codep("MCRS1", "DMAP1")
p <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(
  p,
  filename = "../Outputs/mcrs1_dmap1.pdf", height = 5, width = 10)

p1 <- plot_gene("BRD8", "DMAP1")
p2 <- plot_codep("BRD8", "DMAP1")
p <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(
  p,
  filename = "../Outputs/brd8_dmap1.pdf", height = 5, width = 10)
