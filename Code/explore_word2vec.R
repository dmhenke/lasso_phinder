# Load R libs ####
library(doc2vec)
library(fgsea)


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

write.table(
  aframe,
  file = "/Users/lukas/Downloads/bioplanet_genesets.txt",
  quote = F, row.names = F, col.names = F, sep = "\t")


# Create doc2vec model ####
text <- txt_clean_word2vec(
  aframe$terms, 
  ascii = TRUE, alpha = TRUE, tolower = TRUE, trim = TRUE)
db <- data.frame(doc_id = aframe$gene, text = text)

model <- paragraph2vec(
  x = db,
  type = "PV-DBOW", dim = 50, iter = 20, 
  min_count = 3, lr = 0.05, threads = 4)


predict(
  model,
  newdata = "MYC",
  type = "embedding", which = "docs")

predict(
  model,
  newdata = "splicing",
  type = "embedding", which = "words")

predict(
  model,
  newdata = c("cancer"),
  type = "embedding", which = "words")

sentences <- c("heat shock protein")
sentences <- setNames(sentences, sentences)
sentences <- strsplit(sentences, split = " ")

similar <- predict(
  model,
  newdata = sentences,
  type = "nearest", which = "sent2doc", top_n = 100)
head(similar[[1]], 10)

aframe$terms[aframe$gene == "BCS1L"]


embedding <- as.matrix(model, which = "docs")

coord <- predict(
  model,
  newdata = "metabolism",
  type = "embedding", which = "words")

embedding <- rbind(embedding, coord)
rownames(embedding)[nrow(embedding)] <- "INPUT"

library(umap)

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


newdoc <- doc2vec(model, "i like busses with a toilet")
word2vec_similarity(emb, newdoc)
