# Shane Ridoux
# 250206
# Information Gain

rm(list=ls())
cat("\014")

library(infotheo)
library(data.table)
# source("/Users/shane/School/CU-Denver/Masters-Project/masters-project/entropy.R")

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

# synergy <- function(data = data, pheno = pheno, snps = c(snp1,snp2)){
#   # H <- H(data, var = pheno, conditional = FALSE)
#   snp1 <- snps[1]
#   snp2 <- snps[2]
#   info.abc <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snps)
#   info.ac <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snp1)
#   info.bc <- H(data, pheno) - H(data, pheno, conditional = TRUE, conditioned.on = snp2)
#   
#   synergy <- info.abc-sum(info.ac, info.bc)
#   return(synergy)
# }
# 
# # synergy(data = data, pheno = "Phenotype_discretized", snps = c("chr19:6718376:G:C", "chr16:31265490:G:A"))
# 
# # 
# H(X, "PHENOTYPE") - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps) - sum(H(X, "PHENOTYPE") - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps[1]),
# H(X, "PHENOTYPE") - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps[2]))
# # mutinformation(X) - condinformation(X[2],X[3], S= X[1])
# entropy(X["PHENOTYPE"]) - condentropy(X["PHENOTYPE"], Y = X[colnames(X) != "PHENOTYPE"]) - sum(entropy(X["PHENOTYPE"]) - condentropy(X["PHENOTYPE"], Y = X[2]),entropy(X["PHENOTYPE"]) - condentropy(X["PHENOTYPE"], Y = X[3]))

# synergy <- function(X = X, pheno = pheno, snps = c(snp1,snp2)){
#   H_D <- entropy(X[pheno])
#   syn <- H_D - condentropy(X[pheno], Y = X[colnames(X) != pheno]) - sum(H_D - condentropy(X[pheno], Y = X[snps[1]]),H_D - condentropy(X[pheno], Y = X[snps[2]]))
#   return(syn)
# }
# 
# library(microbenchmark)
# result <- microbenchmark(
#   func1 = synergy_fast(X, "PHENOTYPE", snps = snps),
#   func2 = synergy(X, "PHENOTYPE", snps = snps),
#   times = 10  # Run each function 10 times for better accuracy
# )
# 
# print(result) 
# 
# result <- microbenchmark(
#   func1 = entropy(X["PHENOTYPE"]),
#   func2 = H(X, "PHENOTYPE"),
#   times = 10  # Run each function 10 times for better accuracy
# )
# 
# synergy_faster <- function(X = X, pheno = pheno, snps = c(snp1,snp2)){
#   H_D <- H(X, pheno)
#   info.abc <- H_D - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps)
#   info.ac <- H_D - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps[1])
#   info.bc <- H_D - H(X, "PHENOTYPE", conditional = T, conditioned.on = snps[2])
#   syn <- info.abc - sum(info.ac,info.bc)
#   return(syn)
# }
# 

# H_D <- entropy(X["PHENOTYPE"],method = "emp")
synergy <- function(X = X, pheno = pheno, snps = c(snp1,snp2), entropy = entropy){
  info.abc <- entropy - condentropy(X[,pheno], Y = X[, snps])
  info.ac <- entropy - condentropy(X[,pheno], Y = X[,snps[1]])
  info.bc <- entropy - condentropy(X[,pheno], Y = X[,snps[2]])
  syn <- info.abc - sum(info.ac,info.bc)
  return(syn)
}

# result <- microbenchmark(
#   func1 = synergy(X, "PHENOTYPE", snps = snps),
#   func2 = synergy_faster(X, "PHENOTYPE", snps = snps, entropy = H_D),
#   times = 100
# )
# # 
# print(result)
# 
# library(peakRAM)
# peakRAM(synergy(X, "PHENOTYPE", snps = snps))
# peakRAM(synergy_fast(X, "PHENOTYPE", snps = snps))
# peakRAM(synergy_faster(X, "PHENOTYPE", snps = snps, entropy = H_D))
