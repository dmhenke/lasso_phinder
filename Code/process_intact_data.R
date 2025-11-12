# Load R libs ####
library(biomaRt)
library(data.table)


# Load UniProt IntAct mutations table ####
mut <- fread("/Users/lukas/Downloads/mutations.tsv")
human <- mut[mut$`Affected protein organism` == "9606 - Homo sapiens"]


# Extract interaction pairs ####
get_ids <- function(intx){
  hits_start <- gregexpr("uniprotkb", intx)[[1]]
  unlist(lapply(hits_start, function(x)
    substr(intx, x + 10, x + 15)))  
}

tmp <- human$`Interaction participants`
intxs <- lapply(tmp, get_ids)


# Create adjacency matrix ####
genes <- unique(unlist(intxs))

matr <- matrix(0, length(genes), length(genes))
rownames(matr) <- colnames(matr) <- genes

lapply(intxs, function(x)
  matr[x, x] <<- matr[x, x] + 1)


# Convert UniProt IDs to gene symbols ####
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genelist <- getBM(
  attributes = c("uniprot_gn_id", "hgnc_symbol"),
  filters = "uniprot_gn_id",
  values = genes,
  mart = mart
)
syms <- genelist$hgnc_symbol[match(rownames(matr), genelist$uniprot_gn_id)]
bad <- which(is.na(syms))

final <- matr[-bad, -bad]
rownames(final) <- colnames(final) <- syms[-bad]


# Save output ####
tmp <- reshape2::melt(final)
tmp <- tmp[tmp$Var1 != tmp$Var2, ]
tmp <- tmp[tmp$value > 0, ]
colnames(tmp) <- c("gene1", "gene2", "combined_score")

intact_network <- tmp

save(
  intact_network,
  file = '/Users/lukas/OneDrive/Documents/GitHub/lasso_phinder/Data/intact_network.RData')


