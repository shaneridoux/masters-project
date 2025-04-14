# SER
# 250331
# Inter-gene SNP Synergy Laplacians

rm(list=ls())
cat("\014")

library(data.table)
library(tidyverse)
library(Matrix)
setwd("/Users/shane/School/CU-Denver/Masters-Project/masters-project")
source("information-gain.R")

source("textme.R")

source("within-gene-functions.R")
# catch variables
# args <- commandArgs(trailingOnly = TRUE)

# Get arguments
# chr <- as.numeric(args[1]) 
chr = 4

# read in genotype data
genotype <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19-annotated-pheno.tsv") %>%
  select(-FID,-IID,-MAT,-PAT,-SEX) %>%
  as.data.frame()

pheno <- genotype["PHENOTYPE"]
H_D <- entropy(pheno, method = "emp")

path <- paste0("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/intrxn-pip_chr",chr,".tsv")
intrxns <- fread(path) %>%
  arrange(desc(PIP))

path2 <- paste0("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/main-pip_chr",chr,".tsv")
main <- fread(path2) %>%
  arrange(desc(pip))

main <- main[which(main$pip>0),]

# read in anno 
anno <- fread("/Users/shane/School/CU-Denver/Masters-Project/anno_file.tsv")

for(k in 1:nrow(intrxns)){
gene1 <- intrxns$gene1[k]
gene2 <- intrxns$gene2[k]
pip <- round(intrxns$PIP[k], 3)

gene1_L <- fread(paste0("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians/", gene1, "_L_matrix.csv")) %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

gene2_L <- fread(paste0("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians/", gene2, "_L_matrix.csv")) %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

L_combined <- bdiag(gene1_L, gene2_L)
L <- as.matrix(L_combined)
rownames(L) <- c(rownames(gene1_L), rownames(gene2_L))
colnames(L) <- c(colnames(gene1_L), colnames(gene2_L))
diag(L) = 0
# now i need off-diagonal blocks 
A_off <- L[1:dim(gene1_L)[1],(dim(gene1_L)[1]+1):dim(L)[1]]

genotype_sub <- genotype %>%
  select(all_of(colnames(L)))

X <- cbind(pheno, genotype_sub)

for (i in seq_len(nrow(A_off))) {
  for (j in seq_len(ncol(A_off))) {
    snp1 <- rownames(A_off)[i]
    snp2 <- colnames(A_off)[j]
    A_off[i, j] <- synergy(X = X, pheno = "PHENOTYPE", snps = c(snp1, snp2), entropy = H_D)
  }
}

L[1:dim(gene1_L)[1],(dim(gene1_L)[1]+1):dim(L)[1]] <- -A_off
L[(dim(gene1_L)[1]+1):dim(L)[1], 1:dim(gene1_L)[1]] <- t(-A_off)
diag(L) <- rowSums(-L)

# write out new laplacian
write.table(L,
            file = paste0("/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/Laplacians/",
                          "chr",chr,"_",gene1, "_", gene2, "_pair_L_matrix.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

# get_network(paste0(gene1,"_", gene2), L, output_dir = "/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res")
cat(paste0("Laplacian constructed for chr",chr,": ",gene1," and ", gene2, "!\n"))


pair <- paste0(gene1," and ", gene2)
A = -L
diag(A) = 0

# Set node size by degree centrality
signs <- sign(A)
weights <- abs(A)
threshold <- quantile(weights, 0.99)
G <- graph_from_adjacency_matrix(weights, mode = "undirected", weighted = TRUE,
                                 diag = FALSE)
# Optional: set edge width proportional to magnitude
E(G)$width <- E(G)$weight * 10  # scale as needed

G_top <- delete_edges(G, E(G)[abs(weight) < threshold])
node_size <- degree(G_top, mode = "all")

# Get node names from the graph
node_names <- V(G_top)$name

# Match exonic functions to node names
exonic_function_vector <- anno$exonic_func[match(node_names, anno$topmed)]

# Assign exonic function as a node attribute
V(G_top)$exonic_function <- exonic_function_vector

V(G_top)$gene <- anno$gene[match(node_names, anno$topmed)]
V(G_top)$node_size <- node_size
V(G_top)$label <- ifelse(V(G_top)$node_size > 0, V(G_top)$name, NA)
snp_pos <- stringr::str_extract(V(G_top)$name, "(?<=:)[0-9]+") %>% as.numeric()
snp_rank <- rank(snp_pos, ties.method = "first")

original_names <- V(G_top)$name
new_names <- seq_along(original_names)
name_dict <- data.frame(original = original_names,
                        renamed = new_names,
                        stringsAsFactors = FALSE)
V(G_top)$name <- new_names
V(G_top)$label <- ifelse(V(G_top)$node_size > 0, V(G_top)$name, NA)

E(G_top)$edge_color <- factor(
  ifelse(E(G_top)$weight * signs[as_edgelist(G_top)] < 0, "-", "+"),
  levels = c("+", "-")  # include both
)

# Create a layout
# graph_layout <- layout_with_fr(G_top)  # Force-directed layout
graph_layout <- cbind(
  x = snp_rank,
  y = jitter(degree(G_top), amount = .5)  # adjust `amount` as needed
)
# Plot the network
g <- ggraph(G_top, layout = graph_layout) +
  geom_edge_link(aes(edge_alpha = weight, color = edge_color), show.legend = c(edge_color = TRUE, edge_alpha = FALSE)) +
  scale_edge_color_manual(
    values = c("-" = "#D73027", "+" = "darkgreen"),
    name = "Epistatic Effect",  # Optional custom title
    labels = c("-" = "Antagonistic", "+" = "Synergistic"),
    drop = FALSE
  ) +
  geom_node_point(aes(size = node_size, color = gene), show.legend = c(size = FALSE, color = TRUE)) +  # Color must be inside aes()
  geom_node_text(aes(label = label), color = "black", repel = TRUE, size = 2.5) +
  theme_void() +
  ggtitle(paste0(pair, " Synergy Network (PIP ~ ",pip,")")) +
  scale_color_manual(values = setNames(c("black", "darkblue"), c(gene1, gene2)),
                     name = "Genes")# Save Graph as Image
cat(paste0("Graph constructed for chr",chr,": ",gene1," and ", gene2, "!\n"))
output_file <- file.path("/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/", paste0("Graphs/chr",chr,"_",gene1,"_",gene2, "_network.png"))
ggsave(output_file, plot = g, width = 8, height = 6, dpi = 300, bg = "white")
write.table(name_dict, file = paste0("/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/Graphs/chr",chr,"_",gene1, "_", gene2, "_SNP_dictionary.tsv"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}

for(k in 1:nrow(main)){
  gene <- main$genes[k]
  pip <- round(main$pip[k], 3)
  
  gene_L <- fread(paste0("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians/", gene, "_L_matrix.csv")) %>%
    as.data.frame() %>%
    column_to_rownames("V1") %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  
  A = -gene_L
  diag(A) = 0
  
  # Set node size by degree centrality
  signs <- sign(A)
  weights <- abs(A)
  threshold <- quantile(weights, 0.99)
  G <- graph_from_adjacency_matrix(weights, mode = "undirected", weighted = TRUE,
                                   diag = FALSE)
  # Optional: set edge width proportional to magnitude
  E(G)$width <- E(G)$weight * 10  # scale as needed
  
  G_top <- delete_edges(G, E(G)[abs(weight) < threshold])
  node_size <- degree(G_top, mode = "all")
  
  # Get node names from the graph
  node_names <- V(G_top)$name
  
  # Match exonic functions to node names
  exonic_function_vector <- anno$exonic_func[match(node_names, anno$topmed)]
  
  # Assign exonic function as a node attribute
  V(G_top)$exonic_function <- exonic_function_vector
  
  V(G_top)$gene <- anno$gene[match(node_names, anno$topmed)]
  V(G_top)$node_size <- node_size
  V(G_top)$label <- ifelse(V(G_top)$node_size > 0, V(G_top)$name, NA)
  original_names <- V(G_top)$name
  snp_pos <- stringr::str_extract(V(G_top)$name, "(?<=:)[0-9]+") %>% as.numeric()
  snp_rank <- rank(snp_pos, ties.method = "first")
  
  new_names <- seq_along(original_names)
  name_dict <- data.frame(original = original_names,
                          renamed = new_names,
                          stringsAsFactors = FALSE)
  V(G_top)$name <- new_names
  V(G_top)$label <- ifelse(V(G_top)$node_size > 0, V(G_top)$name, NA)
  
  E(G_top)$edge_color <- factor(
    ifelse(E(G_top)$weight * signs[as_edgelist(G_top)] < 0, "-", "+"),
    levels = c("+", "-")  # include both
  )
  
  # Create a layout
  # graph_layout <- layout_with_fr(G_top)  # Force-directed layout
  # graph_layout <- cbind(
  #   x = snp_pos,
  #   y = jitter(rep(0, length(snp_pos)), amount = 1)  # or use degree, etc.
  # )
  graph_layout <- cbind(
    x = snp_rank,
    y = jitter(degree(G_top), amount = .5)  # adjust `amount` as needed
  )
  
  # Plot the network
  g <- ggraph(G_top, layout = graph_layout) +
    geom_edge_link(aes(edge_alpha = weight, color = edge_color), show.legend = c(edge_color = TRUE, edge_alpha = FALSE)) +
    scale_edge_color_manual(
      values = c("-" = "#D73027", "+" = "darkgreen"),
      name = "Epistatic Effect",  # Optional custom title
      labels = c("-" = "Antagonistic", "+" = "Synergistic"),
      drop = FALSE
    ) +
    geom_node_point(aes(size = node_size, color = exonic_function), show.legend = c(size = FALSE, color = TRUE)) +  # Color must be inside aes()
    geom_node_text(aes(label = label), color = "black", repel = TRUE, size = 2.5) +
    theme_void() +
    ggtitle(paste0(gene, " Synergy Network (PIP ~ ",pip,")")) +
    scale_color_manual(values = c("nonsynonymous SNV" = "darkred", 
                                  "synonymous SNV" = "steelblue", 
                                  "startloss" = "purple",
                                  "stopgain" = "red",
                                  "stoploss" = "darkgreen",
                                  "unknown" = "gray",
                                  "."= "black"))# Save Graph as Image
  cat(paste0("Graph constructed for chr",chr,": ",gene, "!\n"))
  output_file <- file.path("/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/", paste0("Graphs/chr",chr,"_",gene, "_network.png"))
  ggsave(output_file, plot = g, width = 8, height = 6, dpi = 300, bg = "white")
  write.table(name_dict, file = paste0("/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/Graphs/chr",chr,"_",gene, "_SNP_dictionary.tsv"),
              sep = "\t",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
}
cat(paste0("------------------ CHR",chr," DONE! ------------------\n"))
