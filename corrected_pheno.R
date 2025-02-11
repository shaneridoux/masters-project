# Shane Ridoux
# 250203
# Phenotype Correction

rm(list=ls())
cat("\014")

library(data.table)
library(dplyr)

path <- "/Users/shane/School/CU-Denver/Masters-Project/plink2/variant_qc/Step3"

gen <- fread(paste0(path,"/final_qc.pvar"))
sam <- fread(paste0(path,"/final_qc.psam"))

# ancestry pcas
path <- "/Users/shane/School/CU-Denver/Masters-Project/plink2/sample_qc/Step5"

pcs <- fread(paste0(path,"/pca_results.eigenvec"))
pcs <- pcs[which(pcs$IID %in% sam$IID),]   

data <- merge(sam,pcs)
data$PHENO <- data$PHENO-1

# Adjust phenotype using linear regression (for residualization)
pc_model <- glm(PHENO ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, 
                family = binomial(link = "logit"), 
                data = data)

# Extract residuals to get population corrected pheno
data$Phenotype_corrected <- residuals(pc_model, type = "response")
hist(data$Phenotype_corrected, breaks = 30, main = "Histogram of Corrected Phenotype")

# discretize pheno
library(infotheo)
data$Phenotype_discretized <- discretize(data$Phenotype_corrected, 
                                         disc = "equalwidth",
                                         nbins = 2)
data$Phenotype_discretized <- data$Phenotype_discretized -1
data$PHENO - data$Phenotype_discretized

# write out data file with corrected pheno
new_data <- data %>%
  select(-PC1, -PC2, -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10)


write.table(new_data, 
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t",
            file = "/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt")

