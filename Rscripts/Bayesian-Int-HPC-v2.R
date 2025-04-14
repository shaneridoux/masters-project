# Shane Ridoux
# 250327
# Bayesian NLint HPC v2

rm(list=ls())
cat("\014")

library(tidyverse)
library(data.table)
library(NLinteraction)

setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
source("textme.R")

args <- commandArgs(trailingOnly = TRUE)
k <- as.numeric(args[1]) # num arrays

# load api
api <- read.table("api.txt")

dir_path <- "gene-summaries/kPCA"
files <- list.files(dir_path, pattern = "gene_summary_.*\\.tsv$", full.names = TRUE)

X <- files %>%
  map(~ fread(.x) %>%
        column_to_rownames(var = "V1") %>%
        mutate(across(everything(), as.numeric))) %>%  # Force numeric columns
  reduce(bind_cols) %>%
  as.matrix()

dim(X)

Y <- fread("residualized-pheno.txt") %>% 
  select(c(2,5,17)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

colnames(Y) <- c("SEX","PHENOTYPE")


cl <- makeCluster(k)  # Use k clusters
registerDoParallel(cl)

results <- foreach(ns = 1:k, .packages = "NLinteraction", .errorhandling = "pass") %dopar% {
  tryCatch({
    NLint(
      Y = Y[, "PHENOTYPE"],
      X = X,
      C = NULL,
      nIter = 20000,
      nBurn = 2,
      thin = 5,
      nChains = 2,
      ns = ns
    )
  }, error = function(e) e)
}

stopCluster(cl)

saveRDS(results, "BayesianInt-Res/raw_results.rds")
saveRDS(NLmod, paste0("BayesianInt-Res/NLmod_ns", best_ns, ".rds"))

# Select best
valid_results <- results[sapply(results, function(x) !inherits(x, "error"))]
waic <- sapply(valid_results, function(mod) mod$waic)
best_ns <- which.min(waic)
NLmod <- valid_results[[best_ns]]
################# Posterior inclusion probabilities

pip = NLmod$MainPIP
genes = as.vector(colnames(X))
gname = data.frame(genes)
post_inc_props = cbind(gname, pip)

write.table(post_inc_props,
            file = "BayesianInt-Res/genes/main-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

#We now look at the matrix of two-way interaction probabilities.

intMat = NLmod$InteractionPIP
colnames(intMat) <- colnames(X)
rownames(intMat) <- colnames(X)
intrxns <- which(intMat > 0, arr.ind = TRUE)
intrxn_df <- data.frame("gene1"=character(),"gene2"=character(),"PIP"=numeric())
for (i in seq_len(nrow(intrxns))) {
  gene1 <- rownames(intMat)[intrxns[i, 1]]
  gene2 <- colnames(intMat)[intrxns[i, 2]]
  value <- intMat[intrxns[i, 1], intrxns[i, 2]]
  
  # Append new row
  intrxn_df <- rbind(intrxn_df, data.frame(gene1 = gene1, gene2 = gene2, PIP = value))
}

write.table(intrxn_df,
            file = "BayesianInt-Res/gene_pairs/intrxn-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

textme(api = api$V1, project = "masters", channel = "bayesian-int", event = "bayesian-int", description = paste0("Bayesian Interaction in complete!"))

