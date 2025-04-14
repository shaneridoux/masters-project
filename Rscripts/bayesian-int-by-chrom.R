# SER
# 250330
# Bayesian NLint genes in chrom singles for HPC

rm(list=ls())
cat("\014")

library(tidyverse)
library(data.table)
library(NLinteraction)

setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
source("textme.R")

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

# Map task_id (1–22) to (chr, seed)
chromosomes <- 1:22
seeds <- 1
grid <- expand.grid(c = chromosomes, s = seeds)
c <- grid[task_id, "c"]
s <- grid[task_id, "s"]

# load api
api <- read.table("api.txt")

dir_path <- "gene_summaries/kPCA"
files <- list.files(dir_path, pattern = "gene_summary_.*\\.tsv$", full.names = TRUE)

X <- files %>%
  map(~ fread(.x) %>%
        column_to_rownames(var = "V1") %>%
        mutate(across(everything(), as.numeric))) %>%  # Force numeric columns
  reduce(bind_cols) %>%
  as.matrix()

Y <- fread("residualized-pheno.txt") %>% 
  select(c(2,5,17)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

colnames(Y) <- c("SEX","PHENOTYPE")

chr_main <- fread("BayesianInt-Res/chr-pairs/main-pip.tsv") %>%
  as.data.frame()

chr_intrxn <- fread("BayesianInt-Res/chr-pairs/intrxn-pip.tsv") %>%
  as.data.frame() %>%
  arrange(desc(PIP)) 

chr_pairs <- chr_main %>%
  select(chrs) 

#anno 
anno <- fread("anno_file_reduced.tsv") %>% 
  as.data.frame() 

#subset feature space of genes by chromosmal pairs
chrs <- as.character(unlist(chr_pairs[c,]))
chrs_num <- gsub("chr", "", chrs)
genes <- anno$gene[anno$chr %in% chrs_num]
X_sub <- X[,genes]

cat("starting NLmod\n")

waic_log_file <- paste0("BayesianInt-Res/NLmod_chr_", c, "_WAIC_log.txt")


  cat(paste0("Running NLmod with seed ", s, "\n"))

  mod <- NLint(Y = Y[,"PHENOTYPE"], X = X_sub, C = NULL,
               nIter=100000, nBurn=50000, thin=5, nChains=2, ns=1)
  cat("WAIC for seed", s, ":", mod$waic, "\n")
  saveRDS(mod, file = paste0("NLmod100000_chr_",c,"_ns", s, ".rds"))
  
  write(paste(c, s, mod$waic, sep = ","), 
        file = waic_log_file, 
        append = TRUE)

cat("done with NLmod\n")

textme(api = api$V1, project = "masters", channel = "bayesian-int", event = "gene pairs by chrom single", description = paste0("Bayesian Interaction for pair ",chrs[1], " is complete!"))


