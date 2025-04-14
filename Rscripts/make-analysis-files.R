# SER
# 250312
# make analysis file


# clear variables and environment
rm(list=ls())
cat("\014")
setwd("/Users/shane/School/CU-Denver/Masters-Project")

library(data.table)
library(tidyverse)

# load genotype/pheno data
genotype <- fread("genotype-matrix-hg19.raw") %>% as.data.frame()
row.names(genotype) <- genotype$IID

# remove repeated "_Alt" from snp names to follow topmed format 
colnames(genotype) <- sub("_[^_]+$", "", colnames(genotype))

# import annotation file
anno <- fread("anno_file.tsv") %>% as.data.frame()
# anno <- anno[which(anno$exonic_func != "."),] # remove "." (NAs) from exonic function

# double check colnames are in topmed column from anno
length(colnames(genotype)[colnames(genotype) %in% anno$topmed])

cols.keep <- c(colnames(genotype[,1:6]),anno$topmed)

# subset to just annotated snps
genotype <- genotype[,colnames(genotype) %in% cols.keep]

snps <- colnames(genotype)[-c(1:6)] # get snps

gene_snps <- split(anno$topmed, anno$gene)
length(gene_snps)  # Total unique genes
# sapply(gene_snps, length)  # SNP count per gene

# filter to genes with more than one snp
gene_snps_filtered <- gene_snps[sapply(gene_snps, length) > 1]
length(gene_snps_filtered)

str(gene_snps_filtered)

"THRA1/BTR" %in% names(gene_snps_filtered)
names(gene_snps_filtered)[names(gene_snps_filtered) == "THRA1/BTR"] <- "THRA1|BTR"

gene_snps_filtered_sorted <- gene_snps_filtered[order(sapply(gene_snps_filtered, length), decreasing = TRUE)]

selected_cols <- grep(";", names(gene_snps_filtered_sorted), value = TRUE)

gene_snps_filtered_sorted <- gene_snps_filtered_sorted[!names(gene_snps_filtered_sorted) %in% selected_cols]

csmd1 <- genotype[,colnames(genotype) %in% gene_snps_filtered_sorted$CSMD1]

write.table(gene_snps_filtered_sorted$CSMD1, "csmd1_snp_list.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

csdm1_snps.keep <- fread("/Users/shane/School/CU-Denver/Masters-Project/csmd1_pruned_strict.prune.in",
                         header = FALSE)

gene_snps_filtered_sorted$CSMD1 <- csdm1_snps.keep$V1

save(gene_snps_filtered_sorted, file = "gene_snp_map_filtered-v2.RData")

write.table(genotype,
            file = "genotype-matrix-hg19-annotated-pheno.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
