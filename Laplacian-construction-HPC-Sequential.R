# SER
# 250316
# Build and save Laplacian for large genes

rm(list=ls())
cat("\014")

library(tidyverse)
library(infotheo)
library(data.table)

# set wd
setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")

# source handmade functions
source("information-gain.R")

source("textme.R")

source("within-gene-functions.R")


# catch variables
args <- commandArgs(trailingOnly = TRUE)

# Get arguments
chunk_size <- as.numeric(args[1]) # 10,000 genes
# chunk_size <- 1
chunk_num <- as.numeric(args[2]) # a number 1-4
# chunk_num <- 1

# load api
api <- read.table("api.txt")


# load analysis file
genotype <- fread("genotype-matrix-hg19-annotated-pheno.tsv") %>% as.data.frame()
anno <- fread("anno_file.tsv") %>% as.data.frame()

# load gene snp map
load("gene_snp_map_filtered.RData")

# filter gene list by chunk
start <- 1 + (chunk_size * (chunk_num - 1))
stop <- min(chunk_size * chunk_num, length(gene_snps_filtered_sorted))

gene_snp_chunk <- gene_snps_filtered_sorted[start:stop]
gene_name <- names(gene_snps_filtered_sorted)[start:stop]  # Get the gene name
snps_sub <- gene_snps_filtered_sorted[[start:stop]]        # Get the corresponding SNPs

H_D <- entropy(genotype["PHENOTYPE"], method = "emp")

# Generate SNP-SNP combinations
combos <- combn(snps_sub, 2, simplify = FALSE)  # SNP1 - SNP2 pairs
self <- lapply(snps_sub, function(x) c(x, x))   # SNP1 - SNP1 pairs
pairs <- c(combos, self)  # Merge both types

# Compute synergy for each SNP pair
syn_results <- lapply(pairs, function(snp_pair) {
  syn_value <- synergy(X = genotype, pheno = "PHENOTYPE", snps = snp_pair, entropy = H_D)
  
  data.frame(
    Gene = gene_name,
    SNP1 = snp_pair[1],
    SNP2 = snp_pair[2],
    Synergy = syn_value,
    stringsAsFactors = FALSE
  )
})

# Combine results into a single dataframe
results <- bind_rows(syn_results)

L <- get_L(gene_name, results, output_dir = "within-gene-syn-res/Laplacians")
get_network(gene_name, L, output_dir = "within-gene-syn-res/Networks")

textme(api = api$V1, project = "masters", channel = "laplacian",
       event = "Laplacian Construction", description = paste0("Genes ", start, "-", stop))
