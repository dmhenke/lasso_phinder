# Load R libs ####
library(biomaRt)
library(data.table)
library(ggplot2)
library(glmnet)


# Load String data ####

ppi <- fread("../Data/9606.protein.links.v12.0.txt")

ppi$protein1 <- gsub("9606.", "", ppi$protein1, fixed = T)
ppi$protein2 <- gsub("9606.", "", ppi$protein2, fixed = T)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol"), mart = mart)

ppi$gene1 <- gene_info$hgnc_symbol[match(ppi$protein1, gene_info$ensembl_peptide_id)]
ppi$gene2 <- gene_info$hgnc_symbol[match(ppi$protein2, gene_info$ensembl_peptide_id)]


# Load DepMap data ####
# load("C:/Users/Dafydd/Documents/Projects/lasso_PPI/global.RData")
load(paste0("../Data/global.RData"))


# Load functions ####
source("./allfunctions.R")


# Define X and y ####
ok_cells <- intersect(
  rownames(rnaseq),
  rownames(kronos))

X <- rnaseq[ok_cells, ]
y <- kronos[ok_cells,]
expressed <- apply(X, 2, function(x) mean(x > 0))
X <- X[, expressed > 0.95]
X <- apply(X, 2, function(x)
  (x - mean(x))/sd(x))


# test phi, corr vs rmse
gene_sample <- c("MYC",sample(colnames(kronos),100))
gene_sample2 <- c("KDM5D","SOX10","FAM50A","RPP25L","PAX8","KRTAP4-11",
                  "EBF1","IRF4","H2BC15","MYB","MDM2","OR4P4"    ,
                  "KRAS","NRAS","HNF1B","OR4C11","EIF1AX","POU2AF1",  
                  "TP63","BRAF","TTC7A","OR4S2")
gene_sample <-gene_sample[gene_sample%in%colnames(kronos)]
gene_sample2 <-gene_sample2[gene_sample2%in%colnames(kronos)]


library("parallel")

phiout <- mclapply(gene_sample2,function(gene){
  save_results_file <- paste0(gene,'.RData')
  if(save_results_file%in%list.files(paste0("../Outputs/"))){
    print(paste('alredy processed',gene))
    NULL
  }else{
    print(paste('Processing',gene))
    y <- y[, gene]
    scores <- get_scores(gene, ppi)
    out <- run_reg_lasso(
      X, y, scores,
      n_folds = 10, phi_range = seq(0, 1, length = 30))
    save(out,file=paste0("../Outputs/",gene,'.RData'))
    return(out)
  }
},mc.cores = 10)
