# Shane Ridoux
# 250205
# Entropy

# rm(list=ls())
# cat("\014")
# library(infotheo)
# library(data.table)

# data <- fread("/Users/shane/School/CU-Denver/Masters-Project/corrected-pheno.txt")
# 
# D <- data$Phenotype_discretized
# D
# 
# # entropy using infotheo
# HofD <- entropy(D, method = "emp")
# HofD
# # conditional entropy
# condentropy(data$Phenotype_discretized, data$SEX, method = "emp")
# 
# #entropy
# pofd1 <- sum(data$Phenotype_discretized == 1)/nrow(data)
# pofd2 <- 1-pofd1
# entropy <- pofd1*log(1/pofd1) + pofd2*log(1/pofd2)
# 
# # conditional entropy D | sex
# p.d1 <- sum(data$Phenotype_discretized==1)/nrow(data)
# p.d0 <- 1-p.d1
# p.s1 <- sum(data$SEX==1)/nrow(data)
# p.s2 <- 1-p.s1
# p.d0gs1 <- table(data$Phenotype_discretized, data$SEX)[1,1]/sum(table(data$Phenotype_discretized, data$SEX)[,1])
# p.d0gs2 <- table(data$Phenotype_discretized, data$SEX)[1,2]/sum(table(data$Phenotype_discretized, data$SEX)[,2])
# p.d1gs1 <- 1 - p.d0gs1
# p.d1gs2 <- 1 - p.d0gs2
# 
# p.d0s1 <- p.d0gs1*p.s1
# p.d0s2 <- p.d0gs2*p.s2
# p.d1s1 <- p.d1gs1*p.s1
# p.d1s2 <- p.d1gs2*p.s2
# 
# entropy2 <- p.d0s1*log(1/p.d0gs1) + p.d0s2*log(1/p.d0gs2) + p.d1s1*log(1/p.d1gs1) + p.d1s2*log(1/p.d1gs2)



# H <- function(data = data, var = var, conditional = c("TRUE","FALSE",NULL), conditioned.on = NULL){
#   conditional <- match.arg(conditional)
#   
#   if(conditional == "FALSE" | is.null(conditional)){
#     p1 <- sum(data[[var]] == 1)/length(data[[var]])
#     p2 <- 1-p1
#     entropy <- p1*log(1/p1) + p2*log(1/p2)
#     return(entropy)
#   }
#   else if(conditional == "TRUE"){
#     # Create contingency table
#     table_ds <- table(data[[var]], data[[conditioned.on]])
#     
#     # Compute marginal probabilities P(X)
#     p_x <- colSums(table_ds) / sum(table_ds)
#     
#     # Compute conditional probabilities P(D | X)
#     p_d_given_x <- sweep(table_ds, 2, colSums(table_ds), FUN = "/")
#     
#     # Compute joint probabilities P(D, X) = P(D | X) * P(X)
#     p_dx <- sweep(p_d_given_x, 2, p_x, FUN = "*")
#     
#     # Compute entropy sum
#     entropy <- sum(
#       ifelse(p_dx > 0, p_dx * log(1 / p_d_given_x), 0), na.rm = TRUE
#     )
#     return(entropy)
#   }
#   return(NULL)
# }
# 
# H(data, "Phenotype_discretized", conditional = "TRUE", conditioned.on = "SEX")



# H <- function(data, var, conditional = FALSE, conditioned.on = NULL) {
#   # Check if the target variable exists
#   if (!(var %in% colnames(data))) {
#     stop(paste("Error: Column", var, "not found in dataset"))
#   }
#   
#   # If not conditional, compute marginal entropy
#   if (!conditional) {
#     p <- table(data[[var]]) / length(data[[var]])
#     entropy <- -sum(p * log(p), na.rm = TRUE)
#     return(entropy)
#   }
#   
#   # If conditional, check if conditioning variables are provided
#   if (is.null(conditioned.on) || length(conditioned.on) == 0) {
#     stop("Error: 'conditioned.on' must be specified when 'conditional' is TRUE")
#   }
#   
#   # Check if all conditioning variables exist
#   missing_vars <- setdiff(conditioned.on, colnames(data))
#   if (length(missing_vars) > 0) {
#     stop(paste("Error: Columns", paste(missing_vars, collapse = ", "), "not found in dataset"))
#   }
#   
#   # Create a contingency table over the target and conditioned variables
#   table_ds <- table(data[[var]], data[conditioned.on])
#   
#   # Convert the table into a data frame for easier handling
#   table_ds <- as.data.frame(table_ds)
#   colnames(table_ds) <- c(var, conditioned.on, "Freq")
#   
#   # Compute joint probabilities P(D, X)
#   table_ds$P_joint <- table_ds$Freq / sum(table_ds$Freq)
#   
#   # Compute conditional probabilities P(D | X)
#   table_ds$P_X <- ave(table_ds$P_joint, table_ds[, conditioned.on], FUN = sum)
#   table_ds$P_D_given_X <- table_ds$P_joint / table_ds$P_X
#   
#   # Compute conditional entropy H(D | X)
#   entropy <- -sum(ifelse(table_ds$P_D_given_X > 0, 
#                          table_ds$P_joint * log(table_ds$P_D_given_X), 0), 
#                   na.rm = TRUE)
#   
#   return(entropy)
# }

H <- function(data, var, conditional = FALSE, conditioned.on = NULL) {
  # Check if the target variable exists
  if (!(var %in% colnames(data))) {
    stop(paste("Error: Column", var, "not found in dataset"))
  }
  
  # If not conditional, compute marginal entropy
  if (!conditional) {
    p <- table(data[[var]]) / length(data[[var]])
    entropy <- -sum(p * log(p), na.rm = TRUE)
    return(entropy)
  }
  
  # If conditional, check if conditioning variables are provided
  if (is.null(conditioned.on) || length(conditioned.on) == 0) {
    stop("Error: 'conditioned.on' must be specified when 'conditional' is TRUE")
  }
  
  # Check if all conditioning variables exist
  missing_vars <- setdiff(conditioned.on, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste("Error: Columns", paste(missing_vars, collapse = ", "), "not found in dataset"))
  }
  
  # Combine the conditioned variables into a single factor
  conditioned_data <- interaction(data[conditioned.on], drop = TRUE)
  
  # Create a contingency table over the target and the combined conditioning variable
  table_ds <- table(data[[var]], conditioned_data)
  
  # Convert the table into a data frame for easier handling
  table_ds <- as.data.frame(table_ds)
  colnames(table_ds) <- c(var, "conditioned", "Freq")
  
  # Compute joint probabilities P(D, X)
  table_ds$P_joint <- table_ds$Freq / sum(table_ds$Freq)
  
  # Compute conditional probabilities P(D | X)
  table_ds$P_X <- ave(table_ds$P_joint, table_ds$conditioned, FUN = sum)
  table_ds$P_D_given_X <- table_ds$P_joint / table_ds$P_X
  
  # Compute conditional entropy H(D | X)
  entropy <- -sum(ifelse(table_ds$P_D_given_X > 0, 
                         table_ds$P_joint * log(table_ds$P_D_given_X), 0), 
                  na.rm = TRUE)
  
  return(entropy)
}

# H(data, var = "Phenotype_discretized")

# H(data, var = "Phenotype_discretized", conditional = TRUE, conditioned.on = c("chr19:6718376:G:C", "chr16:31265490:G:A"))
# condentropy(data$Phenotype_discretized, Y = data[,c("chr19:6718376:G:C","chr16:31265490:G:A")])
