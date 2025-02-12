# Shane Ridoux
# 250210
# Within Gene Synergy

rm(list=ls())
cat("\014")
set.seed(12)

library(infotheo)
library(data.table)
library(igraph)
library(stringi)
library(dplyr)
library(ggplot2)
library(ggraph)
# BiocManager::install("minet")
library(minet)
source("/Users/shane/School/CU-Denver/Masters-Project/masters-project/entropy.R")
source("/Users/shane/School/CU-Denver/Masters-Project/masters-project/information-gain.R")

data <- fread("/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt")
load("/Users/shane/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/data/complement/daisy/input/analysis_file_3levels.RData")

data <- merge(data, analysis_file, by.x = "IID", by.y = "ID")
data <- as.data.frame(data) # gets rid of data.table class

anno <- as.data.frame(fread("/Users/shane/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/data/complement/daisy/results/Step4-results-all-endpoints-adj-sex-10y.txt"))

snps <- colnames(data)[21:27] # get snps
# rm <- c("chr1:207131555:T:C","chr6:31946247:T:A")
# snps <- snps[-which(snps %in% rm)]
combos <- combn(snps, 2) # make combos of snps to test
self <- matrix(sort(rep(snps,2)),2) # make snpA - snpA pair to test
pairs <- matrix(c(combos,self),2) # combine combos with self rep snp pairs 
syn=c()
SNP1=c()
SNP2=c()
# compute bivariate synergy between snps given pheno
for(i in 1:dim(pairs)[2]){
  snps <- pairs[,i]
  bi <- synergy(data = data, pheno = "Phenotype_discretized", snps = snps)
  SNP1 = c(SNP1, snps[1])
  SNP2 = c(SNP2, snps[2])
  syn = c(syn, bi)
}
bisyn=data.frame(SNP1,SNP2,syn)

## Make a square matrix from the dataframe of bivariate synergy
# get names for row and columns
nameVals <- sort(unique(unlist(bisyn[1:2])))
# construct 0 matrix of correct dimensions with row and column names
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
# fill in the matrix with matrix indexing on row and column names
myMat[as.matrix(bisyn[c("SNP1", "SNP2")])] <- bisyn[["syn"]]
myMat[as.matrix(bisyn[c("SNP2", "SNP1")])] <- bisyn[["syn"]]
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

gene<-density<-mean_dist<-transitivity<-edge_dens<-vertex_con<-edge_con<-NULL
  
############################################################
## select meaningfull edges using minet pkg and summarise graph
###########################################################

#--------------------------------------------------------------------------
# mrnet: Maximum Relevance Minimum Redundancy
graph_mrnet = mrnet(LD)
  
datgraph_mrnet = graph_from_adjacency_matrix(graph_mrnet, mode = "undirected", weighted = TRUE,
                                               diag = FALSE)
#remove loops
datgraph_mrnet = simplify(datgraph_mrnet, remove.multiple=TRUE, remove.loops=TRUE)
  
# gene[[file]] = file
# -----------------------------   summaries of interest
# density = edge_density(datgraph_mrnet,loop=FALSE)  #Density
  
# mean_dist = mean_distance(datgraph_mrnet)  #Average Path Length
  
# transitivity = transitivity(datgraph_mrnet)    #Clustering Coefficeint
  
# edge_dens = edge_density(datgraph_mrnet, loops=F) #number of edges/no.of posible edges
  
# vertex_con = vertex_connectivity(datgraph_mrnet) #number of edges/no.of posible edges
# edge_con = edge_connectivity(datgraph_mrnet) #number of edges/no.of posible edges

# snpbetw_centr = betweenness(datgraph_mrnet, directed=F, weights=NA)
# snpbetw_centr = data.frame(snpbetw_centr)
##dim(snpbetw_centr)
# snpsbetw_centrDF <- tibble::rownames_to_column(snpbetw_centr, "SNP")


# graphxx = data.frame(gene,density,mean_dist,transitivity,edge_dens,vertex_con,edge_con)
# graphxx = data.frame(density,mean_dist,transitivity,edge_dens,vertex_con,edge_con)
  


# graphxx = setNames(graphxx, c("gene", "density", "mean_dist", "transitivity", "edge_dens", "vertex_con", "edge_con"))
# graphxx = setNames(graphxx, c("density", "mean_dist", "transitivity", "edge_dens", "vertex_con", "edge_con"))
  
  
# Set node size by degree centrality
node_size <- degree(datgraph_mrnet, mode = "all")
  
# Create a layout
graph_layout <- layout_with_fr(datgraph_mrnet)  # Force-directed layout

# Example: Creating a mapping between SNP names and labels
node_labels <- c("chr19:6718376:G:C" = "C3_R102G", 
                 "chr1:207131555:T:C" = "C4BPA_I300T",
                 "chr1:207473117:G:A" = "CD21_S639N",
                 "chr6:31946403:G:A" = "CFB_R32Q",
                 "chr1:196673103:G:A" = "CFH_V62I",
                 "chr6:31946247:T:A" = "CFB_L9H",
                 "chr16:31265490:G:A" = "ITGAM_R77H")

# Replace node names with labels in the igraph object
V(datgraph_mrnet)$label <- node_labels[V(datgraph_mrnet)$name]

# Plot the network
ggraph(datgraph_mrnet, layout = graph_layout) +
  geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
  geom_node_point(aes(size = node_size), color = "steelblue") +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Complement SNP Synergy Network")
  