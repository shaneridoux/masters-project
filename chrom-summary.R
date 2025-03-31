# SER
# 250329
# Summarize Chrom

rm(list=ls())
cat("\014")

library(tidyverse)
library(infotheo)
library(data.table)
# library(doParallel)
library(kernlab)


# set wd
# setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
setwd("/Users/shane/School/CU-Denver/Masters-Project/masters-project")
# source handmade functions
source("information-gain.R")

source("textme.R")

source("within-gene-functions.R")

# load api
api <- read.table("api.txt")

setwd("/Users/shane/School/CU-Denver/Masters-Project")
# load analysis file
dir_path <- "HPC-res/gene-summaries/kPCA"
files <- list.files(dir_path, pattern = "gene_summary_.*\\.tsv$", full.names = TRUE)

X <- files %>%
  map(~ fread(.x) %>%
        column_to_rownames(var = "V1") %>%
        mutate(across(everything(), as.numeric))) %>%  # Force numeric columns
  reduce(bind_cols) %>%
  as.matrix()

anno <- fread("anno_file.tsv") 
anno$gene[which(anno$gene == "THRA1/BTR")] <- "THRA1|BTR"
anno <- anno %>% filter(gene %in% colnames(X)) %>%
  select(gene,chr) %>%
  distinct(gene, .keep_all = TRUE) %>% 
  as.data.frame()

pheno <- fread("/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt") %>%
  select(PHENO) %>%
  as.data.frame()

pc_scores <- vector("list", length = 22)
for (i in 1:22){
  Xtmp <- X[,which(colnames(X) %in% anno$gene[which(anno$chr==i)])]

  X_scaled <- scale(Xtmp) %>% as.data.frame()
  kpca_result <- kpca(~., data = X_scaled, kernel = "rbfdot", features = 3)

  rotated_scores <- rotated(kpca_result)
  eigvals <- kpca_result@eig  # or eigenvalues(kpca_result)
  weights <- eigvals[1:ncol(rotated_scores)] / sum(eigvals[1:ncol(rotated_scores)])
  chr_summary <- as.vector(rotated_scores %*% weights)

  pc_scores[[i]] <- chr_summary
  
  # pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  # pc_scores[[i]] <- pca_result$x[, 1]  # First PC scores
}

Y <- fread("residualized-pheno.txt") %>%
  select(Phenotype_corrected)
  as.data.frame()
chrom_summary <- vector("list", length = 22)
for (i in 1:22) {
  Xtmp <- X[, colnames(X) %in% anno$gene[anno$chr == i]]
  X_scaled <- scale(Xtmp) %>% as.data.frame()
  # Fit Gaussian Process Regression
  gp_model <- gausspr(x = X_scaled, y = Y, kernel = "rbfdot")
  
  # Use fitted values as chromosome-level summary
  chrom_summary[[i]] <- predict(gp_model)
}

str(pc_scores)
chrom_summary <- do.call(cbind, chrom_summary)
colnames(chrom_summary) <- paste0("chr", 1:22)

str(chrom_summary)
write.table(chrom_summary, file = "chromosome_summary.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
