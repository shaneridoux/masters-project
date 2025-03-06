# Shane Ridoux
# 250306
# diffusion kPCA

rm(list=ls())
cat("\014")


library(tidyverse)
# load api
api <- read.table("/Users/shane/School/CU-Denver/Masters-Project/masters-project/api.txt")

# load genotype/pheno data
genotype <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw")
genotype <- as.data.frame(genotype) # gets rid of data.table class
row.names(genotype) <- genotype$IID

# remove repeated "_Alt" from snp names to follow topmed format 
colnames(genotype) <- sub("_[^_]+$", "", colnames(genotype))

# import annotation file
anno <- as.data.frame(fread("/Users/shane/School/CU-Denver/Masters-Project/anno_file.tsv"))
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
# genes <- names(gene_snps_filtered)
genes <- c("A1CF","A4GNT;DBR1")
# for each gene, get laplacian and summarize using diffusion kpca
laplacians <- list()
genotype_gene <- list()
path <- "/Users/shane/School/CU-Denver/Masters-Project/within-gene-syn-res/Laplacians/"
for(gene in genes){
  laplacian <- fread(paste0(path, gene, "_L_matrix.csv")) %>%
    as.data.frame() %>%  # Convert from data.table to data.frame
    column_to_rownames(var = "V1") %>%  # Set first column as row names
    mutate_all(as.numeric) %>%  # Ensure all values are numeric
    as.matrix()  # Convert to matrix
  
  snps <- gene_snps_filtered[[gene]]
  genotype_gene[[gene]] <- genotype[,snps]
  laplacians[[gene]] <- list(laplacian,genotype)
  
}

# define diffusion kernel PCA
pcs = 10 
Beta=seq(0,10,0.1)

mydata <- list()
for(i in 1:length(Beta))
{
  beta <- Beta[[i]]
  KG1<-as.matrix(expm(-beta*laplacians[[gene]]))

  mydata[[i]] <- KG1
  names(mydata)[[i]] <- paste0("KG_beta_", i) #Each kernel should have a unique name as required in mixKernel package
}

X <- mydata
ind_beta <- 1/length(X)
#Method2: Calculate every matrix with the equal weight and sum it up
meta.kernel2 <- lapply(as.list(1:length(X)), function(x){
  X[[x]]*ind_beta })
#Get a linear combination of all the matrix to get the final kernel
beta_avg <- Reduce("+",meta.kernel2) #This is our new diffusion parameter


##************************************************************************
## Construct the diffusion kernels with the average diffusion parameter
##************************************************************************
N<-dim(dataxy)[1] # # nb of individuals               
KG<-as.matrix(beta_avg)
datax<-as.matrix(dataxy[,-1])

trdata <- t(datax)

