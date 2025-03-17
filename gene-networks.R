# Shane Ridoux
# 250310
# gene networks

rm(list=ls())
cat("\014")

source("/Users/shane/School/CU-Denver/Masters-Project/masters-project/information-gain.R")
library(tidyverse)
library(igraph)

main <- fread("/Users/shane/School/CU-Denver/Masters-Project/Bayesian-Interaction-Res/main-pip.tsv") %>%
  as.data.frame()

intrxn <- fread("/Users/shane/School/CU-Denver/Masters-Project/Bayesian-Interaction-Res/intrxn-pip.tsv") %>%
  as.data.frame()


summaries <- fread("/Users/shane/School/CU-Denver/Masters-Project/gene_summaries/gene_summary.tsv") %>%
  column_to_rownames(var="V1") %>%
  as.data.frame()
  
pheno <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw") %>% 
  select(c(2,6)) %>%
  column_to_rownames(var = "IID") %>%
  as.data.frame()

df <- summaries
df$pheno <- pheno$PHENOTYPE

# calc synergy
intrxn$synergy <- NA
H_D <- entropy(df$pheno)
for (i in 1:nrow(intrxn)){
  gene1 <- intrxn$gene1[i]
  gene2 <- intrxn$gene2[i]
  
  df_sub <- df %>% 
    select(all_of(c("pheno",gene1,gene2)))
  
  df_sub_disc <- df_sub
  
  df_sub_disc[[gene1]] <- cut(df_sub[[gene1]], breaks = 359, labels = FALSE)
  df_sub_disc[[gene2]] <- cut(df_sub[[gene2]], breaks = 359, labels = FALSE)
  df_sub_disc[[gene1]] <- as.numeric(cut(jitter(df_sub[[gene1]], amount = 1e-6),
                                         breaks = unique(quantile(df_sub[[gene1]], probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
                                         include.lowest = TRUE))
  
  df_sub_disc[[gene2]] <- as.numeric(cut(jitter(df_sub[[gene2]], amount = 1e-6),
                                         breaks = unique(quantile(df_sub[[gene2]], probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
                                         include.lowest = TRUE))
  intrxn$synergy[i] <- synergy(X = df_sub_disc, pheno = "pheno", snps = c(gene1,gene2), entropy = H_D)
}
entropy(df_sub$pheno)
# Create an edge list for the graph
edges <- as.matrix(intrxn[, c("gene1", "gene2")])
gene_network <- graph_from_edgelist(edges, directed = FALSE)
E(gene_network)$weight <- intrxn$PIP
plot(gene_network, 
     vertex.label.cex = 0.8, 
     edge.width = E(gene_network)$weight * 10,  # Scale edges by PIP
     main = "Gene Interaction Network")

