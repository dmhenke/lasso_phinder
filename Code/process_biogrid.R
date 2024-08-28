# Load R libs ####
library(data.table)
library(expm)
library(Matrix)


# Load BIOGRID data ####
dat <- fread(
  "/Users/lukas/Downloads/BIOGRID-ORGANISM-Homo_sapiens-4.4.234.mitab.txt")


# Extract relevant columns ####
geneA <- dat$`Alt IDs Interactor A`
geneA <- unlist(lapply(geneA, function(x)
  strsplit(x, "|", fixed = T)[[1]][2]))
geneA <- gsub("entrez gene/locuslink:", "", geneA, fixed = T)

geneB <- dat$`Alt IDs Interactor B`
geneB <- unlist(lapply(geneB, function(x)
  strsplit(x, "|", fixed = T)[[1]][2]))
geneB <- gsub("entrez gene/locuslink:", "", geneB, fixed = T)


# Create sparse matrix ####
genes <- unique(c(geneA, geneB))
m <- Matrix(nrow = length(genes), ncol = length(genes), data = 0, sparse = TRUE)
rownames(m) <- colnames(m) <- genes

lapply(1:nrow(dat), function(x)
  m[geneA[x], geneB[x]] <<- 1)

save(m, file = "/Users/lukas/Downloads/biogrid_m.RData")
load("/Users/lukas/Downloads/biogrid_m.RData")

# Exponentiate to find neighbors of neighbors ####
diag(m) <- 1

adj_matrix <- m
result <- adj_matrix

for(k in 2:5) {
  result <- result + (adj_matrix %*% adj_matrix)
  adj_matrix <- adj_matrix %*% adj_matrix
}

row_sums <- rowSums(result)
d <- result / row_sums

mean(d == 0)