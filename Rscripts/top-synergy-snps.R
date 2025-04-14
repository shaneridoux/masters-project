# SER
# 250331
# get top synergy pairs of SNPs for top gene-pairs

rm(list=ls())
cat("\014")

library(data.table)
library(tidyverse)

dir <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians"
files <- list.files(dir, full.names = TRUE)

# names of files to extract (top genes)
intrxn <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/intrxn.tsv")
main <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/main.tsv")

top_intrxn <- intrxn %>%
  filter(PIP > quantile(PIP, 0.99))

top_main <- main %>%
  filter(PIP > quantile(PIP, 0.99))

files_main <- paste0(dir, "/", top_main$gene, "_L_matrix.csv")
files_intrxn <- paste0(dir, "/",top_intrxn$chr,"_", top_intrxn$gene1, "_",top_intrxn$gene2,"_pair_L_matrix.csv")

load("/Users/shane/School/CU-Denver/Masters-Project/gene_snp_map_filtered-v2.RData")

extract_synergies <- function(file, gene, snp_map) {
  L <- tryCatch({
    mat <- read_csv(file, col_names = FALSE, show_col_types = FALSE)
    mat <- as.matrix(type.convert(mat, as.is = TRUE))
    storage.mode(mat) <- "numeric"
    mat
  }, error = function(e) {
    message("Failed to read: ", file)
    return(NULL)
  })
  
  snps <- snp_map[[gene]]
  
  if (is.null(snps)) {
    message("Missing SNPs for ", gene)
    return(NULL)
  }
  
  # Fix off-by-one: remove first row/col if needed
  if (nrow(L) == length(snps) + 1) {
    L <- L[-1, -1]
  }
  
  if (length(snps) != nrow(L)) {
    message("Still mismatch: ", gene, " | SNPs: ", length(snps), ", Matrix rows: ", nrow(L))
    return(NULL)
  }
  
  synergy_df <- which(upper.tri(L), arr.ind = TRUE) %>%
    as_tibble() %>%
    mutate(
      SNP1 = snps[row],
      SNP2 = snps[col],
      Synergy = -L[cbind(row, col)],
      Gene = gene
    ) %>%
    select(Gene, SNP1, SNP2, Synergy)
  
  return(synergy_df)
}

all_synergies <- map2_dfr(files_main, top_main$gene, ~ extract_synergies(.x, .y, gene_snps_filtered_sorted))

extract_pair_synergies <- function(file, gene1, gene2, snp_map) {
  # Attempt to read and coerce the matrix
  L <- fread(file) %>% 
      column_to_rownames("V1") %>%
      as.matrix()

  
  snps1 <- snp_map[[gene1]]
  snps2 <- snp_map[[gene2]]
  
  n1 <- length(snps1)
  n2 <- length(snps2)
  
  # Get interaction block between SNPs of gene1 and gene2
  synergy_df <- expand.grid(
    row = 1:n1,
    col = (n1+1):(n1 + n2)
  ) %>%
    as_tibble() %>%
    mutate(
      SNP1 = snps1[row],
      SNP2 = snps2[col - n1],
      Synergy = -L[cbind(row, col)],
      Gene1 = gene1,
      Gene2 = gene2
    ) %>%
    select(Gene1, Gene2, SNP1, SNP2, Synergy)
  
  return(synergy_df)
}

all_pair_synergies <- pmap_dfr(
  list(files_intrxn, top_intrxn$gene1, top_intrxn$gene2),
  ~ extract_pair_synergies(..1, ..2, ..3, gene_snps_filtered_sorted)
)

top_pairs <- all_pair_synergies %>%
  filter(Synergy > quantile(Synergy, 0.99))
top_pairs$type <- "interaction"

top_main <- all_synergies %>%
  filter(Synergy > quantile(Synergy, 0.99))

top_main$Gene2 <- top_main$Gene
colnames(top_main) <- c("Gene1","SNP1","SNP2","Synergy","Gene2")
top_main <- top_main[,c(1,5,2,3,4)]
top_main$type <- "main"

top <- rbind(top_main, top_pairs)

# top_snp_intrxn <- top %>%
#   filter(type == "interaction") %>%
#   group_by(Gene1, Gene2) %>%
#   filter(Synergy > quantile(Synergy, 0.99))
# 
# top_snp_main <- top %>%
#   filter(type == "main") %>%
#   group_by(Gene1, Gene2) %>%
#   filter(Synergy > quantile(Synergy, 0.99))

top_snps <- top %>%
  group_by(type) %>%
  group_by(Gene1, Gene2) %>%
  filter(Synergy > quantile(Synergy, 0.99))

write.table(top_snps, file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/top-snps/top_snps.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
