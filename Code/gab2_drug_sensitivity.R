# Load R libs ####
library(ggplot2)


# Load data ####
drug <- read.csv(
  "/storage/thinc/projects/resources/depmap/data/drug_sensitivity/depmap_drug_auc.csv",
  row.names = 1)
drug_info <- read.csv(
  "/storage/thinc/projects/resources/depmap/data/drug_sensitivity/depmap_drug_metadata.csv")

load("/storage/thinc/projects/resources/depmap/data/global.RData")


# Associate drug sensitivity with CNV ####
ok <- intersect(rownames(cnv), rownames(drug))

gab2 <- cnv[ok, "GAB2"]

res <- apply(drug[ok,], 2, function(x) try(cor.test(x, gab2)))
names(res) <- colnames(drug)

res <- do.call(rbind, lapply(res, function(x)
  as.numeric(x[c("estimate", "p.value")])))
colnames(res) <- c("coef", "pval")

res <- data.frame(drug_info, res)


# Make barplot ####
tmp <- res[res$data_set == "CTD^2", ]
tmp <- tmp[tmp$coef < 0, ]
tmp <- tmp[sort.list(tmp$pval), ]

tmp$drug_name <- factor(
  tmp$drug_name, levels = tmp$drug_name[order(tmp$pval)])

tmp$egfr <- "no"
tmp$egfr[grep("EGFR", tmp$target)] <- "yes"

ggplot(tmp[1:50, ], aes(coef, drug_name, fill = egfr)) +
  labs(
    x = "Correlation coefficient",
    y = "Top 50 drug associations",
    fill = "EGFR target"
  ) +
  scale_fill_manual(
    values = c("grey", "darkred")
  ) +
  geom_bar(stat = 'identity') +
  theme_classic()

ggsave(
  filename = "/storage/thinc/tmp/gab2_drug_barplot.pdf",
  height = 8, width = 5)


# Make scatter plot ####
subm <- data.frame(
  cnv[ok, c("EGFR", "GAB2")], 
  drug = drug[ok, "Drug.sensitivity.AUC..CTD.2..AFATINIB..CTRP.606135."])

ggplot(subm[subm$EGFR > 0.9 & subm$EGFR < 1.1, ], 
       aes(GAB2 > 1.5, drug)) +
  labs(
    title = "DepMap",
    y = "Afatinib sensitivitiy") +
  geom_jitter(width = 0.1) + 
  geom_boxplot(outlier.colour = NA, alpha = 0.1) + 
  stat_compare_means() + 
  theme_classic()

ggsave(
  filename = "/storage/thinc/tmp/gab2_afatinib_plot.pdf",
  height = 6, width = 4)


drug_info[grep("PRMT5", drug_info$target),]

