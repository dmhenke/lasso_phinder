# Load R libs ####
library(doc2vec)
library(fgsea)
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
  type = "PV-DBOW", dim = 50, iter = 20, 
  min_count = 3, lr = 0.05, threads = 4)


# Play around ####
predict(
  model,
  newdata = "28974",
  type = "embedding", which = "docs")

predict(
  model,
  newdata = "splicing",
  type = "embedding", which = "words")

predict(
  model,
  newdata = c("cancer"),
  type = "embedding", which = "words")

sentences <- c("heat shock response")
sentences <- setNames(sentences, sentences)
sentences <- strsplit(sentences, split = " ")

similar <- predict(
  model,
  newdata = sentences,
  type = "nearest", which = "sent2doc", top_n = 100)
similar <- similar[[1]]
similar <- similar[sort.list(-similar$similarity), ]
similar$gene <- gene_list$SYMBOL[match(similar$term2, gene_list$ENTREZID)]
head(similar, 50)

aframe$terms[aframe$gene == "BCS1L"]


# Creat UMAP for a given word ####
gene_list <- read.csv(
  "/Users/lukas/OneDrive/Miko/THINC/projects/cmap/gene_symbols.csv")

embedding <- as.matrix(model, which = "docs")

coord <- predict(
  model,
  newdata = "splicing",
  type = "embedding", which = "words")

embedding <- rbind(embedding, coord)
rownames(embedding)[nrow(embedding)] <- "INPUT"

uData <- umap(embedding)
subm <- data.frame(
  gene = rownames(embedding), uData$layout)

ggplot(subm, aes(x = X1, y = X2)) +
  geom_point() +
  geom_text(
    data = subm[subm$gene == "FANCD2", ],
    aes(label = gene, color = "red")) +
  geom_text(
    data = subm[subm$gene == "INPUT", ],
    aes(label = gene, color = "blue")) +
  theme_minimal()

distances <- proxy::dist(
  uData$layout,
  t(data.matrix(uData$layout[nrow(uData$layout), ]))
)

aframe <- data.frame(
  dist = distances,
  entrez = subm$gene,
  gene = gene_list$SYMBOL[match(subm$gene, gene_list$ENTREZID)]
)
head(aframe[sort.list(aframe$dist), ], 40)

output <- head(aframe$gene[sort.list(aframe$dist)], 50)
write(
  output, 
  file = "/Users/lukas/Downloads/tmp.txt",
  sep = "\n")


