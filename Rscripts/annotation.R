# Shane Ridoux
# 250210
# make analysis file
rm(list=ls())
cat("\014")


library(data.table)

genotype <- as.data.frame(fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19.raw"))
head(genotype[,1:10])

# Extract SNP names (excluding first columns like FID, IID, etc.)
snp_names <- colnames(genotype)[-(1:6)]  

# Convert into a data frame (split by `:`)
snp_df <- data.table(
  SNP = snp_names,
  SNP2 = gsub("chr", "", snp_names),  # Remove "chr" if present
  CHR = sub(":.*", "", snp_names),  # Extract chromosome (before the first ":")
  POS = sub("^[^:]+:([^:]+):.*", "\\1", snp_names)  # Extract numeric position
)

# Remove "chr" prefix to match Ensembl format
snp_df$SNP2 <- sub("_.*", "", snp_df$SNP2)

anno <- as.data.frame(fread("/Users/shane/School/CU-Denver/Masters-Project/final_qc_hg19_anno.hg19_multianno.txt"))
anno$topmed <- paste0(anno$Chr,":",anno$Start,":",anno$Ref,":",anno$Alt)


matching_snps <- anno$topmed[anno$topmed %in% snp_df$SNP2]

length(matching_snps)

anno <- anno[which(anno$topmed %in% matching_snps),c(1:12,26)]

anno <- anno[,c(13,1:12)]

table(anno$ExonicFunc.refGene)
# missense <- anno[which(anno$ExonicFunc.refGene == "nonsynonymous SNV"),]
colnames(anno) <- c("topmed","chr","start","end","ref","alt","function","gene",
                        "gene_detail","exonic_func","amino_acid_change","rsid","euro_maf")
anno_file <- data.frame("gene"=anno$gene,
                        "rsid"=anno$rsid,
                        "topmed"=anno$topmed,
                        "chr"=anno$chr,
                        "location"=anno$start,
                        "ref"=anno$ref,
                        "alt"=anno$alt,
                        "euro_maf"=as.numeric(anno$euro_maf),
                        "amino_acid_change"=anno$amino_acid_change,
                        "func"=anno$`function`,
                        "exonic_func"=anno$exonic_func)

write.table(anno_file,
            file = "/Users/shane/School/CU-Denver/Masters-Project/anno_file.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
