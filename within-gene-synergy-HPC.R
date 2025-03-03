# Shane Ridoux
# 250210
# Within Gene Synergy for HPC

# This script makes the within-gene synergy graphs and saves the Laplacians for 
# downstream analysis (kPCA and then bayesian interaction)

# clear variables and environment
rm(list=ls())
cat("\014")

# set seed
set.seed(12)

# install pacakges
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
# BiocManager::install("minet")
library(minet)
setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
# source handmade functions
source("information-gain.R")
source("textme.R")

# catch variables
args <- commandArgs(trailingOnly = TRUE)
# Get arguments
chunk_size <- args[1] # 10,000 genes
chunk_num <- args[2] # a number 1-4

# load api
api <- read.table("/Users/shane/School/CU-Denver/Masters-Project/masters-project/api.txt")

# load genotype/pheno data
genotype <- fread("genotype-matrix-hg19.raw")
genotype <- as.data.frame(genotype) # gets rid of data.table class
row.names(genotype) <- genotype$IID

# remove repeated "_Alt" from snp names to follow topmed format 
colnames(genotype) <- sub("_[^_]+$", "", colnames(genotype))

# import annotation file
anno <- as.data.frame(fread("anno_file.tsv"))
# anno <- anno[which(anno$exonic_func != "."),] # remove "." (NAs) from exonic function

# double check colnames are in topmed column from anno
length(colnames(genotype)[colnames(genotype) %in% anno$topmed])

cols.keep <- c(colnames(genotype[,1:6]),anno$topmed)

# subset to just annotated snps
genotype <- genotype[,colnames(genotype) %in% cols.keep]

snps <- colnames(genotype)[-c(1:6)] # get snps

gene_snps <- split(anno$topmed, anno$gene)
length(gene_snps)  # Total unique genes
sapply(gene_snps, length)  # SNP count per gene

# filter to genes with more than one snp
gene_snps_filtered <- gene_snps[sapply(gene_snps, length) > 1]

results <- list()  # Store results for each gene
H_D <- entropy(genotype["PHENOTYPE"], method = "emp")
log_file <- "syn_progress_log.txt" # make log file for keeping track of progress
for(gene in names(gene_snps_filtered)){
  write(paste(Sys.time(), "- Processing:", gene), file = log_file, append = TRUE)
  snps_sub <- gene_snps_filtered[[gene]]
  combos <- combn(snps_sub, 2) # SNP A - SNP B
  self <- matrix(sort(rep(snps_sub, 2)), 2) # SNP A - SNP A
  pairs <- matrix(c(combos, self), 2)  # Merge both types of pairs
  num_pairs <- ncol(pairs)
  
  # Preallocate vectors for performance (avoid slow concatenation)
  SNP1 <- character(num_pairs)
  SNP2 <- character(num_pairs)
  syn <- numeric(num_pairs)
  
  syn_results <- lapply(seq_len(ncol(pairs)), function(i) {
    snp_pair <- pairs[, i]  # Extract SNP pair
    syn_value <- synergy(X = genotype, pheno = "PHENOTYPE", snps = snp_pair, entropy = H_D)
    
    # Return a data frame with SNP names and synergy value
    data.frame(
      SNP1 = snp_pair[1],
      SNP2 = snp_pair[2],
      Synergy = syn_value,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine results into a single data frame
  results[[gene]] <- do.call(rbind, syn_results)
}
bisyn <- bind_rows(results, .id = "Gene")
textme(api = api$V1,
       project = "masters",
       channel = "within-gene",
       event = "Synergy Calculation",
       description = "Bivariate synergies have been calculated for all 35898 genes!"
)



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
  datgraph_mrnet = simplify(datgraph_mrnet, remove.multiple=TRUE, remove.loops=TRUE)
  
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

# Enable progress handlers
handlers(global = TRUE)  # Ensures progress messages appear
handlers("txtprogressbar")  # Uses a text-based progress bar

# Get the list of genes
gene_list <- names(results)

# Number of cores to use
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE", unset = detectCores() - 1))

# Run parallel processing
network_results <- mclapply(gene_list, function(g) {
  cat("Processing:", g, "\n")  # Keep log messages for tracking progress
  
  gene_network(gene = g, bisyn = bisyn, 
               output_dir = "within-gene-syn-graphs")
}, mc.cores = num_cores)

# Combine graphxx and snpsbetw_centrDF into a single data frame
network_df <- purrr::map_dfr(network_results, function(x) {
  if (is.null(x$graphxx) || is.null(x$snpsbetw_centrDF)) return(NULL)  # Skip if missing
  
  # Merge gene network data with SNP betweenness centrality
  df <- x$snpsbetw_centrDF %>%
    mutate(gene = x$graphxx$gene,
           density = x$graphxx$density,
           mean_dist = x$graphxx$mean_dist,
           transitivity = x$graphxx$transitivity,
           edge_dens = x$graphxx$edge_dens,
           vertex_con = x$graphxx$vertex_con,
           edge_con = x$graphxx$edge_con)
  
  return(df)
})
textme(api = api$V1,
       project = "masters",
       channel = "within-gene",
       event = "Network Building",
       description = "Networks have been created for all 35898 genes!"
)
# Check the structure of the new combined data frame
str(network_df)

# Display the first few rows
head(network_df)

############# Write out Network Summaries and Laplacians ##################
net.summary <- unique(network_df[,-c(1,2)])
str(net.summary)
write.table(net.summary,
            file = "within-gene-syn-res/network_summary.tsv",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


net.betwn <- data.frame("gene"=network_df$gene,
                        "snp"=network_df$SNP,
                        "betwn"=network_df$snpbetw_centr)
str(net.betwn)
write.table(net.betwn,
            file = "within-gene-syn-res/network_betweeness.tsv",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


# Define directory to save L matrices
output_dir <- "within-gene-syn-res/Laplacians"

# Loop through each element in network_results
for (i in seq_along(network_results)) {
  # Extract the L matrix
  L_matrix <- network_results[[i]]$L
  
  # Check if L exists
  if (!is.null(L_matrix)) {
    # Extract gene name
    gene_name <- network_results[[i]]$graphxx$gene[1] 
    
    # Define file path (saving as CSV)
    file_path <- file.path(output_dir, paste0(gene_name, "_L_matrix.csv"))
    
    # Save L matrix as CSV
    write.csv(L_matrix, file_path, row.names = TRUE)
  }
}
textme(api = api$V1,
       project = "masters",
       channel = "within-gene",
       event = "Laplacian Construction",
       description = "Lacplacians have been saved for all 35898 genes!"
)


write.table(bisyn,
            file = "within-gene-syn-res/bisyn.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
