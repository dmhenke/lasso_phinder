# Load R libs ####
library(biomaRt)
library(data.table)
library(ggplot2)
library(glmnet)
library(parallel)

# Set work folder ####
wrkfldr <- "/mnt/data/user/david/lasso/"
# Load data ####
# Load Depmap data ####
load(paste0(wrkfldr,"../Data/global.RData"))

ppi <- fread(paste0(wrkfldr,"9606.protein.links.v12.0.txt"))
ppi$protein1 <- gsub("9606.", "", ppi$protein1, fixed = T)
ppi$protein2 <- gsub("9606.", "", ppi$protein2, fixed = T)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"), mart = mart)

ppi$gene1 <- gene_info$hgnc_symbol[match(ppi$protein1, gene_info$ensembl_peptide_id)]
ppi$gene2 <- gene_info$hgnc_symbol[match(ppi$protein2, gene_info$ensembl_peptide_id)]

ok_cells <- intersect(
  rownames(rnaseq),
  rownames(kronos))

X <- rnaseq[ok_cells, ]
y <- kronos[ok_cells,]
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))
gene_sample <- c("MYC",sample(colnames(kronos),100))
           
res <- mclapply(gene_sample,function(gene){
  y <- y[, gene]
  scores <- get_scores(gene, ppi)
  out <- run_reg_lasso(
    X, y, scores,
    n_folds = 10, phi_range = seq(0, 1, length = 30))
  save(out,file=paste0(wrkfldr,gene,".RData"))
},mc.cores=10)
