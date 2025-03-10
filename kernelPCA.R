# Shane Ridoux
# 250306
# diffusion kPCA

rm(list=ls())
cat("\014")


library(tidyverse)
library(RSpectra)
# load api
api <- read.table("/Users/shane/School/CU-Denver/Masters-Project/masters-project/api.txt")

# load genotype/pheno data
genotype <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw") %>% as.data.frame()

row.names(genotype) <- genotype$IID

# remove repeated "_Alt" from snp names to follow topmed format 
colnames(genotype) <- sub("_[^_]+$", "", colnames(genotype))

# import annotation file
anno <- fread("/Users/shane/School/CU-Denver/Masters-Project/anno_file.tsv") %>% as.data.frame()
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



# names of genes
path <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians/"
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
for (gene in genes) {
  # Load and process Laplacian matrix
  laplacian <- fread(paste0(path, gene, "_L_matrix.csv")) %>%
    as.data.frame() %>%  
    column_to_rownames(var = "V1") %>%
    mutate_all(as.numeric) %>%
    as.matrix()
  
  # Extract SNPs for the gene
  snps <- gene_snps_filtered[[gene]]
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
            "/Users/shane/School/CU-Denver/Masters-Project/gene_summaries/gene_summary.tsv",
            sep = "\t",
            col.names = T,
            row.names = T)
