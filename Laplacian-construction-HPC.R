# SER
# 250316
# Build and save Laplacian for smaller genes in batches

rm(list=ls())
cat("\014")

library(tidyverse)
library(infotheo)
library(data.table)
library(doParallel)


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
# chunk_size <- 100
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

# filter gene list by chunk
if(chunk_size ==1){
  start <- chunk_num
  stop <- chunk_num
} else if (chunk_size == 100){
  start <- 12 + (chunk_num-1)*(chunk_size)
  stop <- min(start + chunk_size -1, 1397)
  num_cores <- 10
} else if (chunk_size == 1000){
  start <- 1398 + (chunk_num-1)*(chunk_size)
  stop <- min(start + chunk_size -1, 11695)
  num_cores <- 10
} else if (chunk_size == 5000){
  start <- 11696 + (chunk_num-1)*(chunk_size)
  stop <- min(start + chunk_size -1, 20859)
  num_cores <- 10
} else {
  stop("Error: Invalid chunk_size. Please check input values.")
}

# for parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)


gene_snp_chunk <- gene_snps_filtered_sorted[start:stop]
gene_names <- names(gene_snp_chunk) 

H_D <- entropy(genotype["PHENOTYPE"], method = "emp")

results <- foreach(gene = gene_names, .packages = c("dplyr", "infotheo"), .combine = bind_rows) %dopar% {
  # Extract SNPs for the current gene
  snps_sub <- gene_snp_chunk[[gene]]
  
  # Generate SNP-SNP combinations
  combos <- combn(snps_sub, 2)  # SNP1 - SNP2 pairs
  self <- matrix(sort(rep(snps_sub, 2)), 2)  # SNP1 - SNP1 pairs
  pairs <- cbind(combos, self)  # Merge both types of pairs
  num_pairs <- ncol(pairs)
  
  # Compute synergy for each SNP pair using a nested foreach loop
  syn_results <- lapply(seq_len(ncol(pairs)), function(i) {
    snp_pair <- pairs[, i]  # Extract SNP pair
    syn_value <- synergy(X = genotype, pheno = "PHENOTYPE", snps = snp_pair, entropy = H_D)
    
    # Return a data frame with SNP names and synergy value
    data.frame(
      Gene = gene,
      SNP1 = snp_pair[1],
      SNP2 = snp_pair[2],
      Synergy = syn_value,
      stringsAsFactors = FALSE
    )
  })
  
  return(syn_results)
}

# Stop the parallel cluster
stopCluster(cl)

# Process Laplacians and Networks
for (gene in gene_names) {
  gene_name <- gene  # Explicitly define gene_name
  
  gene_results <- results[results$Gene == gene, ]
  
  if (nrow(gene_results) > 0) {
    L <- get_L(gene_name, gene_results, output_dir = "within-gene-syn-res/Laplacians")
    get_network(gene_name, L, output_dir = "within-gene-syn-res/Networks")
  }
}


textme(api = api$V1, project = "masters", channel = "laplacian",
       event = "Laplacian Construction", description = paste0("Genes ", start, "-", stop))
