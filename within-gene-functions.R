# SER
# 250316
# functions for getting Laplacian and constructin ggraph
library(minet)
library(igraph)
library(ggraph)
get_L <- function(gene, bisyn, output_dir){
  gene_data <- bisyn[bisyn$Gene == gene, ]
  
  nameVals <- sort(unique(c(gene_data$SNP1, gene_data$SNP2)))
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(gene_data[c("SNP1", "SNP2")])] <- gene_data$Synergy
  myMat[as.matrix(gene_data[c("SNP2", "SNP1")])] <- gene_data$Synergy
  LD<-myMat
  diag(LD) = 0 #make diagonal zero i.e no info between the same snp
  
  
  ## Make the Diffusion Laplacian matrix
  D<-diag(rowSums(LD))
  Laplacian<-as.matrix(D-LD)
  write.csv(Laplacian,
            file = paste0(output_dir,"/",gene_name,"_L_matrix.csv"),
            row.names = T,
  )
  return(Laplacian) 
}

get_L_chr <- function(chr, bisyn, output_dir){
  gene_data <- bisyn[bisyn$Gene == gene, ]
  
  nameVals <- sort(unique(c(gene_data$SNP1, gene_data$SNP2)))
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(gene_data[c("SNP1", "SNP2")])] <- gene_data$Synergy
  myMat[as.matrix(gene_data[c("SNP2", "SNP1")])] <- gene_data$Synergy
  LD<-myMat
  diag(LD) = 0 #make diagonal zero i.e no info between the same snp
  
  
  ## Make the Diffusion Laplacian matrix
  D<-diag(rowSums(LD))
  Laplacian<-as.matrix(D-LD)
  write.csv(Laplacian,
            file = paste0(output_dir,"/",gene_name,"_L_matrix.csv"),
            row.names = T,
  )
  return(Laplacian) 
}

get_network <- function(gene, Laplacian, output_dir){
  LD <- -Laplacian
  diag(LD) = 0
  
  density<-mean_dist<-transitivity<-edge_dens<-vertex_con<-edge_con<-NULL
  
  ## select meaningfull edges using minet pkg and summarise graph
  
  #--------------------------------------------------------------------------
  # mrnet: Maximum Relevance Minimum Redundancy
  graph_mrnet = mrnet(LD)
  # graph_mrnet = LD
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
    ggtitle(paste(gene, "Synergy Network")) +
    scale_color_manual(values = c("nonsynonymous SNV" = "darkred", 
                                  "synonymous SNV" = "steelblue", 
                                  "startloss" = "purple",
                                  "stopgain" = "red",
                                  "stoploss" = "darkgreen",
                                  "unknown" = "gray",
                                  "."= "black"))
  # Save Graph as Image
  output_file <- file.path(output_dir, paste0("Graphs/",gene, "_network.png"))
  ggsave(output_file, plot = g, width = 8, height = 6, dpi = 300, bg = "white")
  
  # Define output file paths
  network_summary_path <- paste0(output_dir, "/summary/", gene, "_network_summary.csv")
  betweenness_path <- paste0(output_dir, "/betweeness/", gene, "_betweenness.csv")
  res <- list(graphxx = graphxx, snpsbetw_centrDF = snpsbetw_centrDF)
  # Write graph summary data
  write.csv(res$graphxx, file = network_summary_path, row.names = TRUE)
  
  # Write betweenness centrality data
  write.csv(res$snpsbetw_centrDF, file = betweenness_path, row.names = TRUE)
  return(res)
}

