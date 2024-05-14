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


# Get gene coordinates ####
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = colnames(cnv),
  mart = mart)

uniq <- names(which(table(gene_info$hgnc_symbol) == 1))


# Define CNV table ####
X_cnv <- cnv
X_cnv <- na.omit(X_cnv)


# Define outcome ####
gene <- "PRMT5"
y <- demeter2[, gene]
y <- y[!is.na(y)]

tissue <- sample_info[names(y), "lineage"]


# Find overlapping cell lines ####
ok_cells <- intersect(names(y), rownames(X_cnv))


# Remove features without variance ####
X_cnv  <- X_cnv[ok_cells, ]
X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]

y <- y[ok_cells]


# Run LASSO ####
scores <- get_scores(gene, ppi)

results_cnv <- run_reg_lasso(
  X_cnv, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))


# Plot results ####
tmp <- results_cnv$betas
tmp$gene <- rownames(tmp)
tmp$correl <- cor(X_cnv, y)[, 1]

tmp <- tmp[tmp$gene %in% uniq, ]
tmp$start <- gene_info$start_position[match(tmp$gene, gene_info$hgnc_symbol)]
tmp$chr <- gene_info$chromosome_name[match(tmp$gene, gene_info$hgnc_symbol)]

chrs <- as.character(1:22)
tmp <- tmp[tmp$chr %in% chrs, ]
tmp$chr <- factor(tmp$chr, levels = chrs)
tmp <- tmp[order(tmp$chr, tmp$start), ]

tmp$rank <- 1:nrow(tmp)

tmp$label <- tmp$gene
tmp$label[tmp$betas_pen == 0] <- NA
tmp$label[tmp$gene == "SF3B1"] <- "SF3B1"

ggplot(tmp, aes(rank, correl, color = chr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_color_manual(values = rep(c("grey", "black"), 11)) +
  theme_classic()

aframe <- data.frame(y, X_cnv[, c("MTAP", "PRMT5")])
ggplot(aframe, aes(MTAP, y)) +
  geom_point(alpha = 0.5) +
  theme_classic()

p1 <- ggplot(tmp[tmp$chr == 14, ], aes(rank, correl, color = chr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = "black") +
  theme_classic()

p2 <- ggplot(tmp[tmp$chr == 14, ], aes(rank, betas_pen, color = chr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = "black") +
  theme_classic()

tmp2 <- reshape2::melt(
  tmp[tmp$chr == 14, ], measure.vars = c("correl", "betas_pen", "betas"))

ggplot(tmp2, aes(rank, value)) +
  facet_wrap(~variable, ncol = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = "black") +
  theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 1)

ggplot(aframe, aes(PRMT5, y)) +
  labs(
    y = "PRMT dependency",
    x = "PRMT5 copy number"
  ) +
  facet_wrap(~ MTAP > 0.7) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  #stat_cor() +
  theme_classic()
