# Shane Ridoux
# 250210
# make analysis file
rm(list=ls())
cat("\014")


library(data.table)

genotype <- as.data.frame(fread("/Users/shane/School/CU-Denver/Masters-Project/plink2/genotype_matrix.raw"))
head(genotype[,1:10])

# Extract SNP names (excluding first columns like FID, IID, etc.)
snp_names <- colnames(genotype)[-(1:6)]  

# Convert into a data frame (split by `:`)
snp_df <- data.table(
  SNP = snp_names,
  SNP2 = gsub("chr", "", snp_names),
  CHR = sub(":.*", "", snp_names),  # Extract chromosome
  POS = sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names)  # Extract position
)

# Remove "chr" prefix to match Ensembl format
snp_df$SNP2 <- sub("_.*", "", snp_df$SNP2)

fuma <- as.data.frame(fread("/Users/shane/School/CU-Denver/Masters-Project/GCST/FUMA_job581100/snps.txt"))

matching_snps <- fuma$uniqID[fuma$uniqID %in% snp_df$SNP2]
length(matching_snps)

