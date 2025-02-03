# Shane Ridoux
# 250119
# GWAS QC

rm(list=ls())
cat("\014")

################################## SET UP ######################################
# plink GWAS
#!/bin/bash

# Set the working directory
WORKDIR="/Users/shane/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/DAISY/genetics/daisy_ask_genetics/final_plink_report_BDC"
OUTDIR="/Users/shane/School/CU-Denver/Masters-Project/plink2"
PLINKDIR="/Users/shane/School/CU-Denver/Masters-Project/plink2"
# Change to the working directory
# cd $WORKDIR

input_data="daisyask_2" #file for genotype data
keep_file="phenotype_ids.txt"
pheno_file="daisyia.pheno"
phenofile <- read.table(paste0(OUTDIR,"/",pheno_file),
                        header = TRUE)
################################## PRE-FILTERING ###############################
# prefilter out bad quality variants (ie variants with high missing genotype rates)
# $PLINKDIR/plink2 --pfile $input_data --keep $OUTDIR/$keep_file --pheno $OUTDIR/$pheno_file --missing variant-only --out $OUTDIR/pre-filter-missing-var

missingness <- fread(paste0(OUTDIR,"/pre-filter-missing-var.vmiss"))

hist(missingness$F_MISS,
     breaks = 50,
     main = "Genotype Missingness Distribution",
     xlab = "Proportion Missing",
     xlim = c(0,1))

summary(missingness$F_MISS)

################################## SAMPLE QC ###################################
# Step 1
# Sample Call Rate (remove samples with high proportion of missing genotypes based on distribution)
# $PLINKDIR/plink2 --pfile $input_data --keep $OUTDIR/$keep_file --pheno $OUTDIR/$pheno_file --missing sample-only --out $OUTDIR/sample_missingness

sample_missingness <- fread(paste0(OUTDIR,"/sample_missingness.smiss"))

hist(sample_missingness$F_MISS,
     breaks = 50,
     main = "Sample Missingness Distribution",
     xlab = "Proportion Missing",
     xlim = c(0,1))

summary(sample_missingness$F_MISS)

# Step 2
# Check for discrpancies in sex
# The most commonly used QC software (PLINK) will call a sample male if the 
# X chromosome homozygosity rate is more than 0.8; a female call is made if this
# estimate is less than 0.2.

# $PLINKDIR/plink --bfile $bfile --keep $OUTDIR/$keep_file --check-sex --out $OUTDIR/sex_check

# Step 3
# Heterozygosity
# The mean genome-wide heterozygosity of a sample is the fraction or the 
# proportion of non-missing genotypes that are heterozygous in relation to all 
# the genotypes. A reasonable approach is to remove samples that are plus or 
# minus 3 standard deviations (SD) from the mean

# $PLINKDIR/plink2 --pfile $input_data --keep $OUTDIR/$keep_file --het --out $OUTDIR/heterozygosity
het_data <- fread(paste0(OUTDIR,"/heterozygosity.het"))

mean_f <- mean(het_data$F,
               na.rm = TRUE)

sd_f <- sd(het_data$F,
           na.rm = TRUE)

# Define thresholds
lower_bound <- mean_f - 3 * sd_f
upper_bound <- mean_f + 3 * sd_f

# Identify outliers
outliers <- subset(het_data, F < lower_bound | F > upper_bound)

# plot
plot(het_data$F,
     main = "Heterozygosity")
abline(h = lower_bound,
       col = "red")
abline(h = upper_bound,
       col = "red")

remove <- outliers[which(outliers$F < -0.1 | outliers$F > 0.1),]
colnames(remove)[1] <- "FID"
remove <- remove[,c("FID","IID")]

write.table(remove,
            paste0(OUTDIR,"/het-outliers-to-remove.txt"),
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

# $PLINKDIR/plink2 --pfile $input_data --keep $OUTDIR/$keep_file --remove $OUTDIR/het-outliers-to-remove.txt --make-pgen --out $OUTDIR/het_filtered_data
table(phenofile[-which(phenofile$FID %in% remove$FID),"PHENO"])

# Step 4
# Relatedness
# LD pruning and remove complex regions like MHC before determining IBD or kinship
# IBD > 0.9 is used ot identify individuals that are duplicated
# IBD > 0.2 is used to identify individuals that are sexond-degree or closer relative

# LD PRUNING IS NOT RECOMMENDED WITH KING TOOL
# Note that KING kinship coefficients are scaled such that duplicate samples have 
# kinship 0.5, not 1. First-degree relations (parent-child, full siblings) correspond 
# to ~0.25, second-degree relations correspond to ~0.125, etc. It is conventional to 
# use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic 
# twins and duplicate samples, ~0.177 to add first-degree relations, etc

# $PLINKDIR/plink2 --pfile $OUTDIR/het_filtered_data --make-king --out $OUTDIR/king

king_matrix <- as.matrix(read.table(paste0(OUTDIR,"/king.king"), header = FALSE))

# $PLINKDIR/plink2 --pfile $OUTDIR/het_filtered_data --king-cutoff 0.177 --out $OUTDIR/king_filtered

# $PLINKDIR/plink2 --pfile $OUTDIR/het_filtered_data --keep $OUTDIR/king_filtered.king.cutoff.in.id --make-pgen --out $OUTDIR/king_filtered

# Step 5
# Ethnicity
# PCA to identify ethnic outliers which can be removed

# Step 6
# Batch Effects
# gross batch effects can be picked up with PCA

# Step 7
# Sequence Specific Checks for Sample Contamination


################################## VARIANT QC ##################################