# Shane Ridoux
# 250310
# Bayesian NLint

rm(list=ls())
cat("\014")

library(tidyverse)
library(data.table)
library(NLinteraction)

dir_path <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/gene-summaries/kPCA"
files <- list.files(dir_path, pattern = "gene_summary_.*\\.tsv$", full.names = TRUE)

X <- files %>%
  map(~ fread(.x) %>%
        column_to_rownames(var = "V1") %>%
        mutate(across(everything(), as.numeric))) %>%  # Force numeric columns
  reduce(bind_cols) %>%
  as.matrix()

dim(X)

# Y <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw") %>% 
#   select(c(2,5,6)) %>%
#   column_to_rownames(var = "IID") %>%
#   as.matrix()

Y <- fread("/Users/shane/School/CU-Denver/Masters-Project/residualized-pheno.txt") %>% 
  select(c(2,5,17)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

colnames(Y) <- c("SEX","PHENOTYPE")

chr_main <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/chr-pairs/main-pip.tsv") %>%
  as.data.frame()

chr_intrxn <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/chr-pairs/intrxn-pip.tsv") %>%
  as.data.frame() %>%
  arrange(desc(PIP)) 

chr_pairs <- chr_intrxn %>%
  select(chr1,chr2) 

#anno 
anno <- fread("/Users/shane/School/CU-Denver/Masters-Project/anno_file_reduced.tsv") %>% 
  as.data.frame() 

#subset feature space of genes by chromosmal pairs
chrs <- as.character(unlist(chr_pairs[1,]))
chrs_num <- gsub("chr", "", chrs)
genes <- anno$gene[anno$chr %in% chrs_num]
X_sub <- X[,genes]

NLmod = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=1)
# NLmod2 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=2)
# NLmod3 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=3)
# NLmod4 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=4)
# NLmod5 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=5)
# NLmod6 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=6)
# NLmod7 = NLint(Y=Y[,"PHENOTYPE"], X=X_sub, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=7)

# waic <- data.frame("ns"=seq(1,7), "waic"=NA)
# (waic$waic <- c(NLmod1$waic,NLmod2$waic,NLmod3$waic,NLmod4$waic,NLmod5$waic,NLmod6$waic,NLmod7$waic))
# (best_ns <- waic$ns[which.min(waic$waic)])

# NLmod <- switch(best_ns,
#                 `1` = NLmod1,
#                 `2` = NLmod2,
#                 `3` = NLmod3,
#                 `4` = NLmod4,
#                 `5` = NLmod5,
#                 `6` = NLmod6,
#                 `7` = NLmod7
# )


################# Posterior inclusion probabilities

pip = NLmod$MainPIP
genes = as.vector(colnames(X))
gname = data.frame(genes)
post_inc_props = cbind(gname, pip)

write.table(post_inc_props,
            file = "/Users/shane/School/CU-Denver/Masters-Project/Bayesian-Interaction-Res/main-pip.tsv",
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
            file = "/Users/shane/School/CU-Denver/Masters-Project/Bayesian-Interaction-Res/intrxn-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
