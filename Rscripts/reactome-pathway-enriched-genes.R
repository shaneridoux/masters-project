# SER
# 250412
# Reactome pathway enriched genes


rm(list=ls())
cat("\014")

library(tidyverse)
library(data.table)

gene_list <- fread("/Users/shane/School/CU-Denver/Masters-Project/top0.01fdr-reactome.csv", header = FALSE) %>%
  unlist(use.names = FALSE) %>%
  unique()

main <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/main-pip_ALL.tsv") %>%
  filter(pip>0) %>%
  arrange(desc(pip)) %>%
  filter(gene %in% gene_list)

intrxn <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/intrxn-pip_ALL.tsv") %>%
  filter(PIP>0) %>%
  arrange(desc(PIP)) %>%
  filter(gene1 %in% gene_list | gene2 %in% gene_list)



# look at dist'n of pips
plot(density(main$pip), main = "Distribution of Main Effect PIPs")
abline(v = quantile(main$pip, .50))

top_main <- main %>%
  filter(pip>quantile(pip, .50))

plot(density(intrxn$PIP), main = "Distribution of Interaction Effect PIPs")
abline(v = quantile(intrxn$PIP, .50))

top_intrxn <- intrxn %>%
  filter(PIP>quantile(PIP, .50))
# take top

top_names <- unique(c(top_main$gene, top_intrxn$gene1, top_intrxn$gene2))
names <- unique(c(main$gene,intrxn$gene1,intrxn$gene2))

load("/Users/shane/School/CU-Denver/Masters-Project/gene_snp_map_filtered-v2.RData")

gene_snps_filtered_sorted_sub <- gene_snps_filtered_sorted[names(gene_snps_filtered_sorted) %in% names]

dir <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians"

files_main <- paste0(dir, "/", top_main$gene, "_L_matrix.csv")
files_intrxn <- paste0(dir, "/chr",top_intrxn$chr,"_", top_intrxn$gene1, "_",top_intrxn$gene2,"_pair_L_matrix.csv")

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

all_synergies <- map2_dfr(files_main, top_main$gene, ~ extract_synergies(.x, .y, gene_snps_filtered_sorted_sub)) %>%
  arrange(desc(Synergy))
head(all_synergies)
plot(density(all_synergies$Synergy))
abline(v = quantile(all_synergies$Synergy, 0.99))

top_syn <- all_synergies %>%
  filter(Synergy > quantile(all_synergies$Synergy, 0.99))

top_main <- main %>%
  filter(gene %in% unique(top_syn$Gene))

top_main <- left_join(top_main, top_syn, by = c("gene" = "Gene"))

top_main$Gene2 <- top_main$gene
colnames(top_main) <- c("Gene1","PIP","CHR", "SNP1","SNP2","Synergy","Gene2")
top_main <- top_main[,c(3,1,7,2,6,4,5)]
top_main$type <- "main"

top_main_first <- top_main[, .SD[1], by = .(Gene1, Gene2)] %>%
  arrange(desc(Synergy))

write.table(top_main_first, file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/top-snps/top_main_snps.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


extract_pair_synergies <- function(file, gene1, gene2, snp_map) {
  # Attempt to read and coerce the matrix
  L <- fread(file) %>% 
    column_to_rownames("V1") %>%
    as.matrix()
  
  chr <-str_extract(file, "(?<=chr)\\d+")
  pip <- top_intrxn$PIP[which(top_intrxn$gene1 == gene1 & top_intrxn$gene2 == gene2)]
    
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
      CHR = chr,
      PIP = round(pip, 3),
      SNP1 = snps1[row],
      SNP2 = snps2[col - n1],
      Synergy = -L[cbind(row, col)],
      Gene1 = gene1,
      Gene2 = gene2
    ) %>%
    select(Gene1, Gene2, CHR, PIP, Synergy, SNP1, SNP2)
  
  return(synergy_df)
}

all_pair_synergies <- pmap_dfr(
  list(files_intrxn, top_intrxn$gene1, top_intrxn$gene2),
  ~ extract_pair_synergies(..1, ..2, ..3, gene_snps_filtered_sorted_sub)
)

top_pairs <- as.data.table(all_pair_synergies)[Synergy > quantile(Synergy, 0.99)]
top_pairs[, type := "interaction"]

top_intrxn_first <- top_pairs[, .SD[1], by = .(Gene1, Gene2)][order(-Synergy)]

top <- rbind(top_main_first, top_pairs)

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
  filter(Synergy > quantile(Synergy, 0.99)) 

setDT(top_snps)  # Convert to data.table if not already

top_snps_first <- top_snps[
  top_snps[, .I[which.max(Synergy)], by = .(Gene1, Gene2)]$V1
]

top_snps_filtered <- top_snps_first %>%
  filter(Gene1 %in% gene_list | Gene2 %in% gene_list)

write.table(top_snps_filtered, file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/top-snps/top_snps.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")



