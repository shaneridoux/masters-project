# Shane Ridoux
# 250210
# Within Gene Synergy for HPC (Ran Sequentially)

# This script makes the within-gene synergy graphs and saves the Laplacians for 
# downstream analysis (kPCA and then bayesian interaction)

# clear variables and environment
rm(list=ls())
cat("\014")

# set seed
set.seed(12)

# install packages
library(tidyverse)
library(infotheo)
library(data.table)
library(igraph)
library(stringi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggraph)
library(parallel)
library(progressr)
library(progress)
library(doParallel)
library(minet)

# set wd
setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")

# source handmade functions
source("information-gain.R")
source("textme.R")

# catch variables
args <- commandArgs(trailingOnly = TRUE)

# Get arguments
chunk_size <- as.numeric(args[1]) # 10,000 genes
# chunk_size <- 1
chunk_num <- as.numeric(args[2]) # a number 1-4
# chunk_num <- 26900

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

# Detect available cores for parallel processing

num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Log file for tracking progress
log_file <- "synergy_progress_log.txt"

results <- list()  # Store results for each gene
H_D <- entropy(genotype["PHENOTYPE"], method = "emp")
# Run parallel processing with foreach
results <- foreach(gene = names(gene_snp_chunk), .packages = c("dplyr", "infotheo"), .combine = bind_rows) %dopar% {
  # Log progress
  write(paste(Sys.time(), "- Processing:", gene), file = log_file, append = TRUE)
  
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
results <- bind_rows(results)

########### function for gene network and Laplacian #####################
gene_network <- function(gene, bisyn, output_dir){
  gene_data <- bisyn[bisyn$Gene == gene, ]
  
  nameVals <- sort(unique(c(gene_data$SNP1, gene_data$SNP2)))
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(gene_data[c("SNP1", "SNP2")])] <- gene_data$Synergy
  myMat[as.matrix(gene_data[c("SNP2", "SNP1")])] <- gene_data$Synergy
  LD<-myMat
  diag(LD) = 0 #make diagonal zero i.e no info between the same snp
  print(LD)
  
  ## Make the Diffusion Laplacian matrix
  #LD<-myMat 
  D<-diag(rowSums(LD))
  print("D")
  print(D)
  Laplacian<-as.matrix(D-LD)
  print("Laplacian ===========================================")
  print(Laplacian)   
  
  density<-mean_dist<-transitivity<-edge_dens<-vertex_con<-edge_con<-NULL
  
  ## select meaningfull edges using minet pkg and summarise graph
  
  #--------------------------------------------------------------------------
  # mrnet: Maximum Relevance Minimum Redundancy
  graph_mrnet = mrnet(LD)
  
  datgraph_mrnet = graph_from_adjacency_matrix(graph_mrnet, mode = "undirected", weighted = TRUE,
                                               diag = FALSE)
  #remove loops
  if (!any(duplicated(as_edgelist(datgraph_mrnet)))) {
    print("No duplicate edges found. Skipping simplify().")
  } else {
    print("Duplicate edges found. Running simplify() without loop removal.")
    datgraph_mrnet <- simplify(datgraph_mrnet, remove.multiple = TRUE, remove.loops = FALSE)
  }
  
  # -----------------------------   summaries of interest
  density = edge_density(datgraph_mrnet,loop=FALSE)  #Density
  
  mean_dist = mean_distance(datgraph_mrnet)  #Average Path Length
  
  transitivity = transitivity(datgraph_mrnet)    #Clustering Coefficeint
  
  edge_dens = edge_density(datgraph_mrnet, loops=F) #number of edges/no.of posible edges
  
  vertex_con = vertex_connectivity(datgraph_mrnet) #number of edges/no.of posible edges
  edge_con = edge_connectivity(datgraph_mrnet) #number of edges/no.of posible edges
  
  snpbetw_centr = betweenness(datgraph_mrnet, directed=F, weights=NA)
  snpbetw_centr = data.frame(snpbetw_centr)
  
  snpsbetw_centrDF <- tibble::rownames_to_column(snpbetw_centr, "SNP")
  
  
  graphxx = data.frame(gene,density,mean_dist,transitivity,edge_dens,vertex_con,edge_con)
  graphxx = setNames(graphxx, c("gene", "density", "mean_dist", "transitivity", "edge_dens", "vertex_con", "edge_con"))
  
  # Set node size by degree centrality
  node_size <- degree(datgraph_mrnet, mode = "all")
  
  # Create a layout
  graph_layout <- layout_with_fr(datgraph_mrnet)  # Force-directed layout
  
  # Get node names from the graph
  node_names <- V(datgraph_mrnet)$name
  
  # Match exonic functions to node names
  exonic_function_vector <- anno$exonic_func[match(node_names, anno$topmed)]
  
  # Assign exonic function as a node attribute
  V(datgraph_mrnet)$exonic_function <- exonic_function_vector
  # Plot the network
  g <- ggraph(datgraph_mrnet, layout = graph_layout) +
    geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
    geom_node_point(aes(size = node_size, color = exonic_function)) +  # Color must be inside aes()
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    ggtitle(paste(gene_data$Gene, "Synergy Network")) +
    scale_color_manual(values = c("nonsynonymous SNV" = "darkred", 
                                  "synonymous SNV" = "steelblue", 
                                  "startloss" = "purple",
                                  "stopgain" = "red",
                                  "stoploss" = "darkgreen",
                                  "unknown" = "gray",
                                  "."= "yellow"))
  # Save Graph as Image
  output_file <- file.path(output_dir, paste0(gene, "_network.png"))
  ggsave(output_file, plot = g, width = 8, height = 6, dpi = 300, bg = "white")
  
  return(list(graphxx = graphxx, snpsbetw_centrDF = snpsbetw_centrDF, L = Laplacian))
}  


#################### Analysis for gene net and laplacian ########################

# Get the list of genes
gene_list <- unique(results$Gene)

# Number of cores to use
num_cores <- 1

# Run parallel processing
network_results <- gene_network(gene = gene_list, bisyn = results, 
               output_dir = "within-gene-syn-graphs")


# Combine graphxx and snpsbetw_centrDF into a single data frame
network_df <- network_results$snpsbetw_centrDF %>%
    mutate(gene = network_results$graphxx$gene,
           density = network_results$graphxx$density,
           mean_dist = network_results$graphxx$mean_dist,
           transitivity = network_results$graphxx$transitivity,
           edge_dens = network_results$graphxx$edge_dens,
           vertex_con = network_results$graphxx$vertex_con,
           edge_con = network_results$graphxx$edge_con)
  

# Check the structure of the new combined data frame
str(network_df)

# Display the first few rows
head(network_df)

############# Write out Network Summaries and Laplacians ##################
net.summary <- unique(network_df[,-c(1,2)])
str(net.summary)
write.table(net.summary,
            file = paste0("within-gene-syn-res/network_summary_",gene_name,".tsv"),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


net.betwn <- data.frame("gene"=network_df$gene,
                        "snp"=network_df$SNP,
                        "betwn"=network_df$snpbetw_centr)
str(net.betwn)
write.table(net.betwn,
            file = paste0("within-gene-syn-res/network_betweeness_",gene_name,".tsv"),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


# Define directory to save L matrices
output_dir <- "within-gene-syn-res/Laplacians"


# Extract the L matrix
L_matrix <- network_results$L


gene_name <- network_results$graphxx$gene[1]
 
# Define file path (saving as CSV)
file_path <- file.path(output_dir, paste0(gene_name, "_L_matrix.csv"))

# Save L matrix as CSV
write.csv(L_matrix, file_path, row.names = TRUE)


write.table(results,
            file = paste0("within-gene-syn-res/results_",gene_name,".tsv"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

textme(api = api$V1,
       project = "masters",
       channel = "within-gene",
       event = "Laplacian Construction",
       description = paste0("Lacplacian has been saved for gene ",start,"!")
)
