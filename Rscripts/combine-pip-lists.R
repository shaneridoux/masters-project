# SER
# 250331
# combine intrxn/main pips into all chrs

rm(list=ls())
cat("\014")

library(data.table)
library(tidyverse)

dir_path <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs"
intrxn_files <- list.files(dir_path, pattern = "intrxn-pip_.*\\.tsv$", full.names = TRUE)
main_files <- list.files(dir_path, pattern = "main-pip_.*\\.tsv$", full.names = TRUE)

intrxn <- intrxn_files %>%
  map(~ fread(.x)) %>%
  map2_dfr(intrxn_files, ~ mutate(.x, chr = gsub(".*chr(\\d+).*", "chr\\1", .y))) %>%
  select(chr, gene1, gene2, PIP) %>%
  arrange(desc(PIP))

main <- main_files %>%
  map(~ fread(.x)) %>%
  map2_dfr(main_files, ~ mutate(.x, chr = gsub(".*chr(\\d+).*", "chr\\1", .y))) %>%
  select(chr, genes, pip) %>%
  arrange(desc(pip))
colnames(main) <- c("chr","gene","PIP")

write.table(main, file = paste0(dir_path, "/main.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(intrxn, file = paste0(dir_path, "/intrxn.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


write.table(main$gene[which(main$PIP>0)], file = paste0(dir_path, "/main-names.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(c(t(intrxn[, c("gene1", "gene2")])), file = paste0(dir_path, "/intrxn-names.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

###################
main <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/main-pip_ALL.tsv") %>%
  filter(pip>0) %>%
  arrange(desc(pip))

intrxn <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/intrxn-pip_ALL.tsv") %>%
  filter(PIP>0) %>%
  arrange(desc(PIP))

# look at dist'n of pips
plot(density(main$pip))
abline(v = quantile(main$pip, .50))

top_main <- main %>%
  filter(pip>quantile(pip, .50))

plot(density(intrxn$PIP))
abline(v = quantile(intrxn$PIP, .50))

top_intrxn <- intrxn %>%
  filter(PIP>quantile(PIP, .50))
# take top

top_names <- unique(c(top_main$gene, top_intrxn$gene1, top_intrxn$gene2))
names <- unique(c(main$gene,intrxn$gene1,intrxn$gene2))

fwrite(data.frame(gene = names), 
       file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/gene_names_ALL.tsv",
       sep = "\t")

fwrite(data.frame(gene = top_names),
       file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/top_gene_names_ALL.tsv",
       sep = "\t")

fwrite(data.frame(gene = unique(c(intrxn$gene1, intrxn$gene2))),
       file = "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/bayesian-intrxn/gene-pairs/gene_names_intrxn.tsv",
       sep = "\t")
