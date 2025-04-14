# SER
# 250401
# Logistic Regression of top snps

rm(list=ls())
cat("\014")

library(data.table)
library(tidyverse)
library(broom)
library(ggraph)
library(igraph)

top_snps <- fread("/Users/shane/School/CU-Denver/Masters-Project/HPC-res/top-snps/top_snps.tsv") %>%
  as.data.frame()

pheno <- fread("/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt")
pheno2 <- read.delim("/Users/shane/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Immunogenetics_T1D/data/kir_hla/daisy/archive/kir_hla_analysis_df.txt")
pheno2 <- pheno2[pheno2$ID %in% pheno$IID,] %>%
  select(ID,dr34,case_ia_ivyo,NHW)
pheno <- left_join(pheno,pheno2, by = c("IID"="ID"))
rm(list=c("pheno2"))
pheno <- pheno %>%
  select(IID, SEX, NHW, dr34, PHENO)

load("/Users/shane/School/CU-Denver/Masters-Project/gene_snp_map_filtered-v2.RData")
# 
# geno <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19-annotated-pheno.tsv") %>%
#   select(-FID,-SEX,-PHENOTYPE,-PAT,-MAT) %>%
#   column_to_rownames("IID")
# 
# # Step 1: Get all unique SNPs from the gene-to-SNP map
# snp_list <- unique(unlist(gene_snps_filtered_sorted))
# 
# # Step 2: Subset geno to only those SNPs (plus any metadata columns like sample ID)
# snp_cols <- intersect(colnames(geno), snp_list)
# geno_subset <- geno %>% select(all_of(snp_cols))
# 
# write.table(geno_subset, file = "/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19-annotated-subset.tsv",
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

geno <- fread("/Users/shane/School/CU-Denver/Masters-Project/genotype-matrix-hg19-annotated-subset.tsv") %>%
  column_to_rownames("V1")

# ancestry PCs
path <- "/Users/shane/School/CU-Denver/Masters-Project/plink2/variant_qc/Step3"

# gen <- fread(paste0(path,"/final_qc.pvar"))
sam <- fread(paste0(path,"/final_qc.psam"))
path <- "/Users/shane/School/CU-Denver/Masters-Project/plink2/sample_qc/Step5"

pcs <- fread(paste0(path,"/pca_results.eigenvec"))
pcs <- pcs[which(pcs$IID %in% sam$IID),]  


library(pwr)

# Define function to convert odds ratio (OR) to Cohen's f²
# Assumes R² ≈ f² / (1 + f²) => f² = R² / (1 - R²)
# Use approximation R² = log(OR)^2 / (π^2 / 3) for logistic regression
or_to_f2 <- function(or) {
  log_or <- log(or)
  r2 <- (log_or^2) / (pi^2 / 3)
  f2 <- r2 / (1 - r2)
  return(f2)
}

# Try for protective OR = 0.67 (i.e. inverse of OR = 1.5)
f2 <- or_to_f2(0.67)

# Set number of predictors
u <- 6  # SNP1, SNP2, SNP1*SNP2, SEX, dr34, PC1
v <- 359 - u - 1  # degrees of freedom

# Desired power
desired_power <- 0.8

# Solve for required alpha (detectable significance level)
power_result <- pwr.f2.test(u = u, v = v, f2 = f2, sig.level = NULL, power = desired_power)
alpha_detectable <- power_result$sig.level

# Use Bonferroni correction: how many tests would allow alpha_detectable to be Bonferroni-adjusted at alpha=0.1
max_tests <- floor(0.1 / alpha_detectable)
max_tests

snps <- top_snps %>%
  arrange(desc(Synergy)) %>%
  head(max_tests)

tidy <- list()

for(i in 1:nrow(snps)){
  pair <- as.character(snps[i,c("SNP1","SNP2")])
  gene_pair <- as.character(snps[i,c("Gene1","Gene2")])
  df <- cbind(pheno, geno[,colnames(geno) %in% pair])
  df <- merge(df, pcs[,c("IID","PC1")], by = "IID")

  formula_str <- paste0("PHENO ~ `", pair[1], "` + `", pair[2], "` + `", 
                      pair[1], "`*`", pair[2], "` + SEX + dr34 + PC1")

  model <- glm(as.formula(formula_str), data = df, family = "binomial")
  
  (tidy_or <- tidy(model, exponentiate = TRUE))
  
  # Add confidence intervals manually (fallback)
  conf <- tryCatch(
    confint.default(model), 
    error = function(e) NULL
  )
  
  # If conf worked, bind CI to tidy results
  if (!is.null(conf)) {
    conf <- exp(conf)  # because you exponentiated the tidy
    tidy_or <- cbind(tidy_or, conf.low = conf[, 1], conf.high = conf[, 2])
  } else {
    tidy_or$conf.low <- NA
    tidy_or$conf.high <- NA
  }
  
  tidy_or <- mutate(tidy_or, SNP1 = pair[1], SNP2 = pair[2], Gene1 = gene_pair[1], Gene2 = gene_pair[2])
  
  tidy[[i]] <- tidy_or
  
}

tidy <- bind_rows(tidy)
tidy <- tidy %>%
  mutate(
    effect_type = case_when(
      str_detect(term, "`:`") ~ "interaction",
      term %in% c("SEX", "dr34", "PC1") ~ "covariate",
      term == "(Intercept)" ~ "intercept",
      TRUE ~ "main"
    )
  )

tidy <- tidy %>%
  filter(effect_type %in% c("interaction","main")) %>%
  arrange(p.value)

# 1. Filter interaction terms
interactions <- tidy %>%
  filter(effect_type == "interaction")

interactions <- interactions %>%
  mutate(p.adj.bh = p.adjust(p.value, method = "BH"))

sig <- tidy %>%
  filter(effect_type == "interaction") %>%
  filter(p.value < 0.1)

# maybe penalized logistic regression
# library(glmnet)
# 
# for (i in 1:nrow(snps)){
#   ints <- paste0("`",snps$SNP1,"`*`",snps$SNP2,"`")
#   mains <- paste0("`",unique(c(snps$SNP1,snps$SNP2)),"`")
#   covs <- "SEX + dr34"
#   formula_str <- paste(c(mains, ints, covs), collapse = " + ")
#   formula <- as.formula(paste("PHENO ~", formula_str))
# }
# 
# # Build the model matrix (automatically handles interactions)
# geno_sub <- geno[,colnames(geno) %in% unique(c(snps$SNP1,snps$SNP2))]
# X <- model.matrix(formula, data = cbind(pheno, geno_sub))[, -1]
# y <- pheno$PHENO
# penalty <- rep(1, ncol(X))
# penalty[which(colnames(X) == "dr34")] <- 0  # Don't penalize dr34
# 
# # Fit logistic regression with elastic net penalty (alpha=1 is LASSO, alpha=0 is Ridge)
# cv_fit <- cv.glmnet(X, y,
#                     family = "binomial",
#                     alpha = 0,
#                     standardize = TRUE,
#                     penalty.factor = penalty)
# 
# # View lambda that minimizes cross-validation error
# best_lambda <- cv_fit$lambda.min
# print(best_lambda)
# 
# # Fit the final model using best lambda
# final_model <- glmnet(X, y, family = "binomial", alpha = 0, lambda = best_lambda)
# 
# # View coefficients
# coef(final_model)


dir <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/within-gene-syn-res/Laplacians"

# Output directory
out_dir <- "/Users/shane/School/CU-Denver/Masters-Project/HPC-res/log-reg"
anno <- fread("/Users/shane/School/CU-Denver/Masters-Project/anno_file.tsv")
for (j in 1:nrow(sig)) {
  row <- sig[j, ]
  
  # Get appropriate Laplacian file
  if (row$Gene1 == row$Gene2) {
    file <- list.files(dir, pattern = paste0(row$Gene1, "_L_matrix.csv"), full.names = TRUE)
  } else {
    file <- list.files(dir, pattern = paste0(".*_", row$Gene1, "_", row$Gene2, "_pair_L_matrix.csv"), full.names = TRUE)
  }
  
  if (length(file) == 0) {
    warning(paste("No Laplacian file found for", row$Gene1, row$Gene2))
    next
  }
  
  # Read Laplacian and prepare adjacency matrix
  L <- fread(file) %>%
    column_to_rownames("V1") %>%
    as.matrix()
  
  A <- -L
  diag(A) <- 0
  signs <- sign(A)
  weights <- abs(A)
  threshold <- quantile(weights, 0.95, na.rm = TRUE)
  
  G <- graph_from_adjacency_matrix(weights, mode = "undirected", weighted = TRUE, diag = FALSE)
  E(G)$width <- E(G)$weight * 10
  G_top <- delete_edges(G, E(G)[abs(weight) < threshold])
  node_size <- degree(G_top)
  
  node_names <- V(G_top)$name
  V(G_top)$exonic_function <- anno$exonic_func[match(node_names, anno$topmed)]
  V(G_top)$gene <- anno$gene[match(node_names, anno$topmed)]
  V(G_top)$node_size <- node_size
  
  # Color edges by sign
  E(G_top)$edge_color <- factor(
    ifelse(signs[as_edgelist(G_top)] < 0, "-", "+"),
    levels = c("+", "-")
  )
  
  # Position layout
  snp_pos <- stringr::str_extract(node_names, "(?<=:)[0-9]+") %>% as.numeric()
  snp_rank <- rank(snp_pos, ties.method = "first")
  graph_layout <- cbind(x = snp_rank, y = jitter(degree(G_top), amount = 0.5))
  
  # Relabel prioritized SNPs
  original_names <- V(G_top)$name
  new_names <- seq_along(original_names)
  name_dict <- data.frame(original = original_names, renamed = new_names, stringsAsFactors = FALSE)
  
  snp1_index <- which(name_dict$original == row$SNP1)
  snp2_index <- which(name_dict$original == row$SNP2)
  
  if (length(snp1_index) == 0 | length(snp2_index) == 0) {
    warning(paste("Missing SNP index for", row$SNP1, row$SNP2))
    next
  }
  
  new_names[snp1_index] <- row$SNP1
  new_names[snp2_index] <- row$SNP2
  V(G_top)$name <- new_names
  V(G_top)$label <- ifelse(V(G_top)$node_size > 0, V(G_top)$name, NA)
  V(G_top)$focal <- V(G_top)$name %in% c(sig$SNP1, sig$SNP2)
  
  # Odds Ratio label coordinates
  coords <- graph_layout
  x1 <- coords[snp1_index, "x"]
  y1 <- coords[snp1_index, "y"]
  x2 <- coords[snp2_index, "x"]
  y2 <- coords[snp2_index, "y"]
  xmid <- (x1 + x2) / 2
  ymid <- (y1 + y2) / 2
  xlab <- xmid + 10
  ylab <- ymid + 2
  
  label_text <- paste0(
    "OR = ", round(row$estimate, 2),
    " [", round(row$conf.low, 2), ", ", round(row$conf.high, 2), "]"
  )
  
  # Plot
  gene1 <- row$Gene1
  gene2 <- row$Gene2
  pip <- unique(round(snps %>% filter(Gene1 == gene1, Gene2 == gene2) %>% pull(PIP), 3))
  if(gene1 == gene2){
    pair <- gene1
  }else{
    pair <- paste0(gene1, " and ", gene2)
  }
  g <- ggraph(G_top, layout = graph_layout) +
    geom_edge_link(aes(edge_alpha = weight, color = edge_color), 
                   show.legend = c(edge_alpha = FALSE, color = TRUE)) +  # Keep only edge color legend
    scale_edge_color_manual(
      values = c("-" = "#D73027", "+" = "darkgreen"),
      name = "Epistatic Effect",
      labels = c("-" = "Negative", "+" = "Positive"),
      drop = FALSE
    ) +
    geom_node_point(aes(size = node_size, color = gene, shape = focal), 
                    stroke = 2, 
                    show.legend = c(size = FALSE, color = TRUE, shape = FALSE)) +  # Keep only node color legend
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 19)) +
    geom_node_text(aes(label = label), 
                   color = "black", 
                   repel = TRUE, 
                   size = 2.5, 
                   box.padding = unit(1.2, "lines"), 
                   point.padding = unit(1.9, "lines")) +
    theme_void() +
    ggtitle(paste0(pair, " Synergy Network (PIP ~ ", pip, ")")) +
    scale_color_manual(values = setNames(c("black", "lightblue"), c(gene1, gene2)), 
                       name = "Genes")
  
  g2 <- g +
    geom_segment(aes(x = xlab, y = ylab, xend = xmid, yend = ymid),
                 arrow = arrow(type = "closed", length = unit(0.2, "cm")),
                 color = "red") +
    annotate("text", x = xlab, y = ylab, label = label_text, size = 3.5, vjust = -0.4, color = "red")
  # Save
  output_file <- file.path(out_dir, paste0(gene1, "_", gene2, "_network.png"))
  ggsave(output_file, plot = g2, width = 8, height = 6, dpi = 300, bg = "white")
  
  write.table(name_dict,
              file = file.path(out_dir, paste0(gene1, "_", gene2, "_SNP_dictionary.tsv")),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  message("✅ Graph saved for: ", gene1, " and ", gene2)
}
