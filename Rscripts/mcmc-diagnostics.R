# SER
# 250404
# Full MCMC diagnostics and PIP export + post-processing

library(data.table)
library(tidyverse)
library(coda)
library(NLinteraction)

rm(list = ls())
cat("\014")

setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
source("textme.R")
# Load API
api <- read.table("api.txt")

# Load gene summaries
dir_path <- "gene_summaries/kPCA"
files <- list.files(dir_path, pattern = "gene_summary_.*\\.tsv$", full.names = TRUE)

X <- files %>%
  map(~ fread(.x) %>%
        column_to_rownames(var = "V1") %>%
        mutate(across(everything(), as.numeric))) %>%
  reduce(bind_cols) %>%
  as.matrix()

# Load phenotype
Y <- fread("residualized-pheno.txt") %>%
  select(c(2, 5, 17)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()
colnames(Y) <- c("SEX", "PHENOTYPE")

# Annotation
anno <- fread("anno_file_reduced.tsv") %>% as.data.frame()

# Create output folders if needed
dir.create("BayesianInt-Res/diagnostics", showWarnings = FALSE, recursive = TRUE)
dir.create("BayesianInt-Res/gene_pairs", showWarnings = FALSE)

diag_file <- "BayesianInt-Res/diagnostics/diagnostics_log.tsv"
if (!file.exists(diag_file)) {
  write.table(data.frame(chr = "chr", psrf = "psrf", ess = "ess"),
              file = diag_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Process each chromosome
for (chr in 1:22) {
  rds_path <- paste0("NLmod_chr_", chr, "_ns1.rds")
  if (!file.exists(rds_path)) {
    cat("Skipping chr", chr, " — no RDS file\n")
    next
  }
  
  genes <- anno$gene[anno$chr == chr]
  if (length(genes) == 0) {
    cat("No genes found for chromosome", chr, "\n")
    next
  }
  
  X_sub <- X[, genes]
  cat("Processing Chromosome", chr, "\n")
  NLmod <- readRDS(rds_path)
  cat("WAIC:", NLmod$waic, "\n")
  
  ## MCMC Diagnostics
  mcmc_betaC <- mcmc.list(
    mcmc(t(NLmod$posterior$betaC[1, , drop = FALSE])),
    mcmc(t(NLmod$posterior$betaC[2, , drop = FALSE]))
  )
  
  # Save traceplot
  png(paste0("BayesianInt-Res/diagnostics/traceplot_chr", chr, ".png"), width = 800, height = 400)
  traceplot(mcmc_betaC, main = paste("Traceplot for betaC, chr", chr))
  dev.off()
  
  # Convergence diagnostics
  diag <- gelman.diag(mcmc_betaC)
  ess <- effectiveSize(mcmc_betaC)
  
  write.table(data.frame(chr = chr, psrf = diag$psrf[1], ess = ess),
              file = diag_file, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  ## Main PIP
  pip <- NLmod$MainPIP
  main_out <- data.frame(gene = colnames(X_sub), pip = pip)
  fwrite(main_out, file = paste0("BayesianInt-Res/gene_pairs/main-pip_", chr, ".tsv"),
         sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ## Interaction PIP
  intMat <- NLmod$InteractionPIP
  colnames(intMat) <- colnames(X_sub)
  rownames(intMat) <- colnames(X_sub)
  
  intrxns <- which(intMat > 0, arr.ind = TRUE)
  intrxn_df <- data.frame(
    gene1 = rownames(intMat)[intrxns[, 1]],
    gene2 = colnames(intMat)[intrxns[, 2]],
    PIP = intMat[intrxns]
  )
  
  fwrite(intrxn_df, file = paste0("BayesianInt-Res/gene_pairs/intrxn-pip_", chr, ".tsv"),
         sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Notify
  textme(api = api$V1, project = "masters", channel = "bayesian-int",
         event = "gene pairs by chrom single",
         description = paste0("Bayesian interaction model for chr", chr, " complete!"))
}

## Merge PIP files after loop
cat("Merging main/intrxn PIPs across chromosomes...\n")

main_all <- rbindlist(lapply(1:22, function(chr) {
  file <- paste0("BayesianInt-Res/gene_pairs/main-pip_", chr, ".tsv")
  if (file.exists(file)) fread(file)[, chr := chr]
}))
fwrite(main_all, "BayesianInt-Res/gene_pairs/main-pip_ALL.tsv", sep = "\t")

intrxn_all <- rbindlist(lapply(1:22, function(chr) {
  file <- paste0("BayesianInt-Res/gene_pairs/intrxn-pip_", chr, ".tsv")
  if (file.exists(file)) fread(file)[, chr := chr]
}))
fwrite(intrxn_all, "BayesianInt-Res/gene_pairs/intrxn-pip_ALL.tsv", sep = "\t")

cat("All done!\n")