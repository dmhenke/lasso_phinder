# Load R libs ####
library(data.table)
library(doc2vec)
library(fgsea)
library(ggplot2)
library(lsa)
library(umap)
library(word2vec)


# Load GO Mol Function ####
paths <- gmtPathways("/Users/lukas/Downloads/GO_Molecular_Function_2023.txt")
names(paths) <- unlist(lapply(names(paths), function(x)
  strsplit(x, "(", fixed = T)[[1]][1]))

allgenes <- unique(unlist(paths))
terms <- unique(names(paths))

matr <- matrix(0, nrow = length(allgenes), ncol = length(terms))
rownames(matr) <- allgenes
colnames(matr) <- terms

lapply(names(paths), function(x)
  matr[paths[[x]], x] <<- 1)

genes <- lapply(rownames(matr), function(x)
  paste(colnames(matr)[matr[x, ] == 1], collapse = ","))
names(genes) <- rownames(matr)


# Save output ####
aframe <- data.frame(
  gene = names(genes),
  terms = unlist(genes))
write.table(
  aframe,
  file = "/Users/lukas/Downloads/go_mf_genesets.txt",
  quote = F, row.names = F, col.names = F, sep = "\t")


# Load BioPlanet ####
tmp <- read.csv("/Users/lukas/Downloads/pathway.csv")

asplit <- split(tmp$PATHWAY_NAME, tmp$GENE_SYMBOL)
asplit <- lapply(asplit, function(x)
  paste(x, collapse = ","))
aframe <- data.frame(
  gene = names(asplit),
  terms = unlist(asplit))


# Save output ####
write.table(
  aframe,
  file = "/Users/lukas/Downloads/bioplanet_genesets.txt",
  quote = F, row.names = F, col.names = F, sep = "\t")


# Load Generif ####
tmp <- fread(
  "/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/generifs_basic.gz", 
  sep = "\t")
tmp <- tmp[tmp[["#Tax ID"]] == 9606, ]
tmp <- tmp[-grep("(HuGE Navigator)", fixed = T, tmp[["GeneRIF text"]]),]
tail(sort(table(tmp[["GeneRIF text"]])))

generif <- tmp
aframe <- data.frame(
  gene = generif[["Gene ID"]],
  terms = generif[["GeneRIF text"]])


# Create doc2vec model ####
text <- txt_clean_word2vec(
  aframe$terms, 
  ascii = TRUE, alpha = TRUE, tolower = TRUE, trim = TRUE)
db <- data.frame(doc_id = aframe$gene, text = text)

model <- paragraph2vec(
  x = db,
  type = "PV-DBOW", dim = 100, iter = 20, 
  min_count = 5, lr = 0.05, threads = 4)


# Load gene annotations ####
gene_list <- read.csv(
  "/Users/lukas/OneDrive/Miko/THINC/projects/cmap/gene_symbols.csv")


# Play around ####
predict(
  model,
  newdata = "1",
  type = "embedding", which = "docs")

predict(
  model,
  newdata = "splicing",
  type = "embedding", which = "words")

predict(
  model,
  newdata = c("cancer"),
  type = "embedding", which = "words")

sentences <- c("breast cancer")
sentences <- setNames(sentences, sentences)
sentences <- strsplit(sentences, split = " ")

similar <- predict(
  model,
  newdata = sentences,
  type = "nearest", which = "sent2doc", top_n = 100)
similar <- similar[[1]]
similar <- similar[sort.list(-similar$similarity), ]
similar$gene <- gene_list$SYMBOL[match(similar$term2, gene_list$ENTREZID)]
head(similar, 30)


# Check some genes ####
gene <- "KLHL22"
gene_list$ENTREZID[gene_list$SYMBOL == gene]

db[db$doc_id == "641384", ]


# Create UMAP for a given word ####
embedding <- as.matrix(model, which = "docs")

input <- "cancer"
coord <- predict(
  model,
  newdata = input,
  type = "embedding", which = "words")

embedding <- rbind(embedding, coord)
rownames(embedding)[nrow(embedding)] <- "INPUT"

uData <- umap(embedding)
subm <- data.frame(
  id = rownames(embedding), uData$layout)

subm$gene <- gene_list$SYMBOL[match(subm$id, gene_list$ENTREZID)]

helicases <- union(
  unique(gene_list$SYMBOL[grep("^DHX", gene_list$SYMBOL)]),
  unique(gene_list$SYMBOL[grep("^DDX", gene_list$SYMBOL)]))

ggplot(subm, aes(x = X1, y = X2)) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  geom_point(alpha = 0.1) +
  geom_point(
    data = subm[subm$gene %in% helicases, ], color = "blue") +
  ggrepel::geom_text_repel(
    data = subm[subm$gene %in% helicases, ],
    aes(label = gene), color = "blue",
    max.overlaps = 100, size = 3) +
  geom_text(
    data = subm[subm$id == "INPUT", ],
    aes(label = input), color = "red", size = 5) +
  theme_classic()


# Calculate cosine distances ####
dist <- apply(embedding, 1 , function(x)
  cosine(x, embedding[nrow(embedding), ]))

subm <- data.frame(
  id = names(dist),
  gene = gene_list$SYMBOL[match(names(dist), gene_list$ENTREZID)],
  dist = dist)


# Source biolasso ####
source("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Code/allfunctions.R")


# Load DepMap data ####
load("/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/global.RData")


# Normalize RNA expression data ####
X_cnv <- cnv
X_cnv <- na.omit(X_cnv)
X_cnv <- X_cnv[, apply(X_cnv, 2, var) > 0]
X <- X_cnv

X <- rnaseq
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))


# Define outcome ####
gene <- "SF3B1"
y <- demeter2[, gene]
y <- y[!is.na(y)]


# Run LASSO ####
ok_cells <- intersect(
  rownames(X),
  names(y)
)

X <- X[ok_cells, ]
y <- y[ok_cells]

scores <- subm$dist
names(scores) <- subm$gene
scores <- (scores - min(scores)) / (max(scores) - min(scores))

results <- run_reg_lasso(
  X, y, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

aframe <- data.frame(
  gene = rownames(results$betas),
  results$betas,
  score = scores[match(rownames(results$betas), names(scores))],
  cor = cor(X, y)[,1])

ggplot(aframe, aes(
  betas, betas_pen, 
  color = cor)) +
  labs(
    x = "Baseline",
    y = "Bio-primed"
  ) +
  ggrepel::geom_text_repel(
    data = aframe[aframe$betas_pen != 0 |
                    aframe$betas != 0, ],
    aes(label = gene, size = score),
    max.overlaps = 20
  ) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red") +
  theme_classic()

subm <- data.frame(
  cnv = X[, "WDR33"],
  dep = y
)
ggplot(subm, aes(x = cnv, y = dep)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


