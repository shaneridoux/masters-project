# Shane Ridoux
# 250306
# diffusion kPCA HPC

rm(list=ls())
cat("\014")


library(tidyverse)
library(RSpectra)
library(data.table)
library(Matrix)
# setwd("/Users/shane/School/CU-Denver/Masters-Project")
setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
source("textme.R")
# catch variables
args <- commandArgs(trailingOnly = TRUE)

# Get arguments
chunk_size <- as.numeric(args[1]) # 10,000 genes
# chunk_size <- 1000
chunk_num <- as.numeric(args[2]) # a number 1-4
# chunk_num <- 4


# load api
api <- read.table("api.txt")

# load analysis file
genotype <- fread("genotype-matrix-hg19-annotated-pheno.tsv") %>% as.data.frame()

# load gene snp map
load("gene_snp_map_filtered-v2.RData")

# names of genes that completed laplacian step
path <- "within-gene-syn-res/Laplacians/"
files <- list.files(path)
genes <- sub("_L_matrix.csv","",files)

# for each gene, get laplacian and summarize using diffusion kpca
laplacians <- list()
genotype_gene <- list()
pcs = 10 
Beta=seq(0,10,0.1)

KL_data <- list()
K0 <- list()
K <- list()
eigs <- list()
gene_summary <- list()

# filter gene list by chunk
start <- 1 + (chunk_size * (chunk_num - 1))
stop <- min(chunk_size * chunk_num, length(genes))

gene_snp_chunk <- genes[start:stop]

for (gene in gene_snp_chunk) {
  # Load and process Laplacian matrix
  print(gene)
  laplacian <- fread(paste0(path, gene, "_L_matrix.csv")) %>%
    as.data.frame() %>%  
    column_to_rownames(var = "V1") %>%
    mutate_all(as.numeric) %>%
    as.matrix()
  
  # Extract SNPs for the gene
  snps <- gene_snps_filtered_sorted[[gene]]
  genotype_gene[[gene]] <- genotype[, snps] %>% as.matrix()
  laplacians[[gene]] <- list(laplacian, genotype_gene[[gene]])
  
  # Diffusion Kernel PCA
  for (i in 1:length(Beta)) {
    beta <- Beta[[i]]
    KL <- as.matrix(expm(beta * laplacians[[gene]][[1]]))
    
    KL_data[[i]] <- KL
    names(KL_data)[[i]] <- paste0("KL_", i)
  }
  
  # Compute average diffusion kernel
  ind_beta <- 1 / length(KL_data)
  meta.kernel2 <- lapply(as.list(1:length(KL_data)), function(x) {
    KL_data[[x]] * ind_beta
  })
  KL_avg <- as.matrix(Reduce("+", meta.kernel2))
  
  ##### Construct K = G * KL_avg * G^T
  K0[[gene]] <- genotype_gene[[gene]] %*% KL_avg %*% t(genotype_gene[[gene]])
  
  # Center K matrix
  K[[gene]] <- scale(K0[[gene]], center = TRUE, scale = FALSE)
  
  # Eigen decomposition (scaled by N)
  N <- dim(K[[gene]])[1]  # Number of individuals
  eigs[[gene]] <- eigs_sym(K[[gene]] / N, k = pcs, which = "LM", sigma = NULL, lower = TRUE, retvec = TRUE)
  
  # Projection: K * eigenvectors
  eigVector <- eigs[[gene]]$vectors
  Yx <- K[[gene]] %*% eigVector
  
  # Store first PC as gene summary
  gene_summary[[gene]] <- Yx[, 1]
}

gene_summary_df <- as.data.frame(do.call(cbind, gene_summary))


# write out summaries
write.table(gene_summary_df,
            paste0("gene_summaries/kPCA/gene_summary_",start,"-",stop,".tsv"),
            sep = "\t",
            col.names = T,
            row.names = T)

textme(api = api$V1, project = "masters", channel = "kernelpca", event = "kernelPCA", description = paste0("KernelPCA on genes ", start,"-",stop," are done!"))
