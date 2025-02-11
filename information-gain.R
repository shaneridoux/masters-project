# Shane Ridoux
# 250206
# Information Gain

rm(list=ls())
cat("\014")

library(infotheo)
library(data.table)
source("/Users/shane/School/CU-Denver/Masters-Project/masters-project/entropy.R")

# data <- fread("/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt")
# load("/Users/shane/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/data/complement/daisy/input/analysis_file_3levels.RData")

# data <- merge(data, analysis_file, by.x = "IID", by.y = "ID")
# data <- as.data.frame(data)
# H(data, var = "Phenotype_discretized", conditional = FALSE)
# # Information Gain
# 
# condentropy(data$Phenotype_discretized, Y = data[,c("chr19:6718376:G:C","chr16:31265490:G:A")])
# H(data, "Phenotype_discretized", conditional = TRUE, conditioned.on = c("chr19:6718376:G:C","chr16:31265490:G:A"))
# 
# 
# info.abc <- H(data, "Phenotype_discretized") - H(data, "Phenotype_discretized", conditional = TRUE, conditioned.on = c("chr19:6718376:G:C","chr16:31265490:G:A"))
# info.ac <- H(data, "Phenotype_discretized") - H(data, "Phenotype_discretized", conditional = TRUE, conditioned.on = "chr16:31265490:G:A")
# info.bc <- H(data, "Phenotype_discretized") - H(data, "Phenotype_discretized", conditional = TRUE, conditioned.on = "chr19:6718376:G:C")
# 
# # synergy <- info.abc-sum(info.ac, info.bc)
# # -interinformation(dataxy[,c(1,7,8)], method = "emp")

synergy <- function(data = data, pheno = pheno, snps = c(snp1,snp2)){
  # H <- H(data, var = pheno, conditional = FALSE)
  snp1 <- snps[1]
  snp2 <- snps[2]
  info.abc <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snps)
  info.ac <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snp1)
  info.bc <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snp2)
  
  synergy <- info.abc-sum(info.ac, info.bc)
  return(synergy)
}

# synergy(data = data, pheno = "Phenotype_discretized", snps = c("chr19:6718376:G:C", "chr16:31265490:G:A"))



