# Shane Ridoux
# 250329
# Bayesian NLint for CHROM

rm(list=ls())
cat("\014")

library(tidyverse)
library(data.table)
library(NLinteraction)

setwd("/scratch/alpine/sridoux@xsede.org/ms-proj")
source("textme.R")

# args <- commandArgs(trailingOnly = TRUE)
# k <- as.numeric(args[1]) # num arrays

# load api
api <- read.table("api.txt")


X <- fread("chromosome_summary.tsv") %>%
  as.matrix()

summary(colMeans(X))

Y <- fread("residualized-pheno.txt") %>% 
  select(c(2,5,17)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

colnames(Y) <- c("SEX","PHENOTYPE")


NLmod1 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=1)
NLmod2 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=2)
NLmod3 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=3)
NLmod4 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=4)
NLmod5 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=5)
NLmod6 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=6)
NLmod7 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=10000, nBurn=2, thin=5, nChains=2, ns=7)

waic <- data.frame("ns"=seq(1,7), "waic"=NA)
waic$waic <- c(NLmod1$waic,NLmod2$waic,NLmod3$waic,NLmod4$waic,NLmod5$waic,NLmod6$waic,NLmod7$waic)
best_ns <- waic$ns[which.min(waic$waic)] 

NLmod <- switch(best_ns,
                `1` = NLmod1,
                `2` = NLmod2,
                `3` = NLmod3,
                `4` = NLmod4,
                `5` = NLmod5,
                `6` = NLmod6,
                `7` = NLmod7
)
NLmod = NLmod1

################# Posterior inclusion probabilities

pip = NLmod$MainPIP
chrs = as.vector(colnames(X))
cname = data.frame(chrs)
post_inc_props = cbind(cname, pip)

write.table(post_inc_props,
            file = "BayesianInt-Res/chr-pairs/main-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

#We now look at the matrix of two-way interaction probabilities.

intMat = NLmod$InteractionPIP
colnames(intMat) <- colnames(X)
rownames(intMat) <- colnames(X)
intrxns <- which(intMat > 0, arr.ind = TRUE)
intrxn_df <- data.frame("chr1"=character(),"chr2"=character(),"PIP"=numeric())
for (i in seq_len(nrow(intrxns))) {
  chr1 <- rownames(intMat)[intrxns[i, 1]]
  chr2 <- colnames(intMat)[intrxns[i, 2]]
  value <- intMat[intrxns[i, 1], intrxns[i, 2]]
  
  # Append new row
  intrxn_df <- rbind(intrxn_df, data.frame(chr1 = chr1, chr2 = chr2, PIP = value))
}

write.table(intrxn_df,
            file = "BayesianInt-Res/chr-pairs/intrxn-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

textme(api = api$V1, project = "masters", channel = "chr",
       event = "chr intrxn", description = "Chromosome Interaction is Complete!")
