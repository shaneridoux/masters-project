# Shane Ridoux
# 250310
# Bayesian NLint

rm(list=ls())
cat("\014")

library(tidyverse)
library(NLinteraction)

X <- fread("/Users/shane/School/CU-Denver/Masters-Project/gene_summaries/gene_summary.tsv") %>% 
  column_to_rownames(var = "V1") %>%
  as.matrix()


Y <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw") %>% 
  select(c(2,5,6)) %>%
  column_to_rownames(var = "IID") %>%
  as.matrix()

NLint <- function (Y = Y, X = X, C = C, nChains = 2, nIter = 10000, nBurn = 2000, 
          thin = 8, c = 0.001, d = 0.001, sigB = "EB", k = 15, ns = 3, 
          alph = 3, gamm = dim(X)[2], threshold = 0.1, intMax = 3, 
          speed = TRUE) 
{
  n = dim(X)[1] # 359
  p = dim(X)[2] # 447
  designC = cbind(rep(1, dim(X)[1]), C) # 359 x 1
  SigmaC = 1000 * diag(dim(designC)[2]) # 1 x 1
  muC = rep(0, dim(designC)[2]) # 1 x 1
  muB = rep(0, ns) # 1 x ns 
  Xstar = array(NA, dim = c(n, p, ns + 1)) # 359 x 447 x (ns+1)
  Xstar[, , 1] = 1 
  for (j in 1:p) {
    Xstar[, j, 2:(ns + 1)] = scale(splines::ns(X[, j], df = ns))
  }
  if (speed == TRUE) {
    if (sigB == "EB") {
      SigMin = NLinteraction:::MCMCmixtureMinSig(Y = Y, X = X, C = C, Xstar = Xstar, 
                                 nPerms = 10, nIter = 500, c = c, d = d, sigBstart = 0.5, 
                                 muB = muB, SigmaC = SigmaC, muC = muC, k = k, 
                                 ns = ns, threshold = threshold)
      print("Finding the empirical bayes estimate of the slab variance")
      SigEstEB = NLinteraction:::MCMCmixtureEB(Y = Y, X = X, C = C, Xstar = Xstar, 
                               nChains = nChains, nIter = nIter, nBurn = nBurn, 
                               thin = thin, c = c, d = d, sigBstart = 0.5, muB = muB, 
                               SigmaC = SigmaC, muC = muC, SigMin = SigMin, 
                               k = k, ns = ns, alph = alph, gamm = gamm, intMax = intMax)
      SigEst = max(SigEstEB, SigMin)
      print("Now fitting the posterior conditional on EB estimate of slab variance")
      posterior = NLinteraction:::MCMCmixture(Y = Y, X = X, C = C, Xstar = Xstar, 
                              nChains = nChains, nIter = nIter, nBurn = nBurn, 
                              thin = thin, c = c, d = d, sigB = SigEst, muB = muB, 
                              SigmaC = SigmaC, muC = muC, k = k, ns = ns, alph = alph, 
                              gamm = gamm, intMax = intMax)
    }
    else {
      print("Fitting the posterior using user selecte value for slab variance")
      posterior = MCMCmixture(Y = Y, X = X, C = C, Xstar = Xstar, 
                              nChains = nChains, nIter = nIter, nBurn = nBurn, 
                              thin = thin, c = c, d = d, sigB = sigB, muB = muB, 
                              SigmaC = SigmaC, muC = muC, k = k, ns = ns, alph = alph, 
                              gamm = gamm, intMax = intMax)
    }
    keep = nBurn + (1:floor((nIter - nBurn)/thin)) * thin
    totalScans = length(keep)
    waic = NLinteraction:::WaicMixture(Y = Y, Xstar = Xstar, designC = designC, 
                       totalScans = totalScans, nChains = nChains, zetaPost = posterior$zeta, 
                       betaList = posterior$beta, betaCPost = posterior$betaC, 
                       sigmaPost = posterior$sigma, n = n, k = k, ns = ns)
    intMean = NLinteraction:::InteractionMatrix(zetaPost = posterior$zeta, 
                                totalScans = totalScans, nChains = 2, p = p, k = k)
    inclusions = NLinteraction:::InclusionVector(zetaPost = posterior$zeta, 
                                 totalScans = totalScans, nChains = 2, p = p, k = k)
  }
  else {
    if (sigB == "EB") {
      SigMin = NLinteraction:::MCMCmixtureMinSig_MH(Y = Y, X = X, C = C, 
                                    Xstar = Xstar, nPerms = 10, nIter = 500, c = c, 
                                    d = d, sigBstart = 0.5, muB = muB, SigmaC = SigmaC, 
                                    muC = muC, k = k, ns = ns, threshold = threshold)
      print("Finding the empirical bayes estimate of the slab variance")
      SigEstEB = NLinteraction:::MCMCmixtureEB_MH(Y = Y, X = X, C = C, 
                                  Xstar = Xstar, nChains = nChains, nIter = nIter, 
                                  nBurn = nBurn, thin = thin, c = c, d = d, sigBstart = 0.5, 
                                  muB = muB, SigmaC = SigmaC, muC = muC, k = k, 
                                  ns = ns, alph = alph, gamm = gamm, SigMin = SigMin, 
                                  probSamp1 = 0.95, intMax = intMax)
      SigEst = max(SigEstEB, SigMin)
      print("Now fitting the posterior conditional on EB estimate of slab variance")
      posterior = MCMCmixture_MH(Y = Y, X = X, C = C, Xstar = Xstar, 
                                 nChains = nChains, nIter = nIter, nBurn = nBurn, 
                                 thin = thin, c = c, d = d, sigB = SigEst, muB = muB, 
                                 SigmaC = SigmaC, muC = muC, intMax = intMax, 
                                 k = k, ns = ns, alph = alph, gamm = gamm, probSamp1 = 0.95)
    }
    else {
      print("Fitting the posterior using user selecte value for slab variance")
      posterior = MCMCmixture_MH(Y = Y, X = X, C = C, Xstar = Xstar, 
                                 nChains = nChains, nIter = nIter, nBurn = nBurn, 
                                 thin = thin, c = c, d = d, sigB = sigB, muB = muB, 
                                 SigmaC = SigmaC, muC = muC, intMax = intMax, 
                                 k = k, ns = ns, alph = alph, gamm = gamm, probSamp1 = 0.95)
    }
    keep = nBurn + (1:floor((nIter - nBurn)/thin)) * thin
    totalScans = length(keep)
    waic = WaicMixture_MH(Y = Y, Xstar = Xstar, designC = designC, 
                          totalScans = totalScans, nChains = nChains, zetaPost = posterior$zeta, 
                          betaList = posterior$beta, betaCPost = posterior$betaC, 
                          sigmaPost = posterior$sigma, n = n, k = k, ns = ns)
    intMean = InteractionMatrix_MH(zetaPost = posterior$zeta, 
                                   totalScans = totalScans, nChains = 2, p = p, k = k)
    inclusions = InclusionVector_MH(zetaPost = posterior$zeta, 
                                    totalScans = totalScans, nChains = 2, p = p, k = k)
  }
  l = list(posterior = posterior, waic = waic, InteractionPIP = intMean, 
           MainPIP = inclusions, ns = ns, k = k, speed = speed)
  return(l)
}

NLmod1 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=1000, nBurn=2, thin=5, nChains=2, ns=1)
# NLmod2 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=2)
# NLmod3 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=3)
# NLmod4 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=4)
# NLmod5 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=5)
# NLmod6 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=6)
# NLmod7 = NLint(Y=Y[,"PHENOTYPE"], X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=7)

# waic <- data.frame("ns"=seq(1,7), "waic"=NA)
# waic$waic <- c(NLmod1$waic,NLmod2$waic,NLmod3$waic,NLmod4$waic,NLmod5$waic,NLmod6$waic,NLmod7$waic)
# best_ns <- waic$ns[which.min(waic$waic)] # 6

NLmod = NLmod1

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
  entropy <- 
  # Append new row
  intrxn_df <- rbind(intrxn_df, data.frame(gene1 = gene1, gene2 = gene2, PIP = value))
}

write.table(intrxn_df,
            file = "/Users/shane/School/CU-Denver/Masters-Project/Bayesian-Interaction-Res/intrxn-pip.tsv",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
