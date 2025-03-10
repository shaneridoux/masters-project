# Shane Ridoux
# 250309
# Bayesian Interaction (Between Gene Synergy)

rm(list=ls())
cat("\014")

library(tidyverse)
library(NLinteraction)

X <- fread("/Users/shane/School/CU-Denver/Masters-Project/gene_summaries/gene_summary.tsv") %>% 
  column_to_rownames(var = "V1") %>%
  as.matrix()


Y <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw") %>% 
  select(c(2,6)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

NLmod2 = NLint(Y=Y, X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=1)

NLmod = NLmod2

##So we can now evaluate the WAIC of each model
print("WAIC for ns=1 model")
print(NLmod2$waic)

################===========================================
#Posterior inclusion probabilities

pdp = NLmod$MainPIP
genes = as.vector(colnames(X))
gname = data.frame(genes)
post_inc_props = cbind(gname, pdp)

## write out gene names in a txt
#writeLines(genenam, "genes_BayesianTrial.txt")

# write.csv(post_inc_props, "post_inc_probs_Iter20k_ns1.csv")
##===========================================================================

####### =======================================================================
#We now look at the matrix of two-way interaction probabilities.

intMat = NLmod$InteractionPIP
# write.csv(intMat, "intMat_2way_iter20kns1.csv")
#### =========================================================================

## =========================================================================
# image
## save plot interactionProbabilities Plot
pdf(file="int_probs_2way_20k_ns1.pdf")
plotInt(NLmod = NLmod)
dev.off()

## ========================================================================

## write out gene names in a txt
writeLines(genenam, "genes_BayesianTrial.txt")

# =================================================================
## save trace plot
pdf(file="Trace_plot_20kns1.pdf")
plot(NLmod2$posterior$sigma[1,], type="l")
lines(NLmod2$posterior$sigma[2,], col=2)
dev.off()

# =================================================================



