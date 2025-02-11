# Shane Ridoux
# kpca 
# 250203

# rm(list=ls())
# cat("\014")
# 
# data("iris")
# head(iris)
# 
# # just take features and make into matrix
# X = as.matrix(iris[,-5])
# 
# # initialize kernel matrix
# N = nrow(X)
# i = N # individuals
# j = N # genes
# K = matrix(data = 0, nrow = i, ncol = j)
# 
# # Compute the RBF kernel: K(i, j) = exp(-gamma * ||xi - xj||^2)
# # for (i in 1:N) {
# #   for (j in 1:N) {
# #     dist_sq <- sum((X[i, ] - X[j, ])^2)
# #     K[i, j] <- exp(-gamma * dist_sq)
# #   }
# # }
# 
# # linear kernel for regular PCA
# K <- X %*% t(X)
# 
# # center kernel
# N <- nrow(K)
# H <- diag(N) - (1/N)*rep(1,N)%*%t(rep(1,N)) # centering matrix
# K_centered <- H%*%K%*%H # removes row means (HK) and col means (KH)
# 
# # eigen decomposition
# eig <- eigen(K_centered)  # Compute eigenvalues and eigenvectors
# eigvals <- eig$values  # Eigenvalues
# eigvecs <- eig$vectors  # Eigenvectors
# 
# # Set negative eigenvalues or very small values (close to zero) to zero
# eigvals[eigvals < 1e-10] <- 0  # Threshold for small eigenvalues (tune as needed)
# 
# # Normalize eigenvectors (same as PCA)
# Z <- as.data.frame(eigvecs %*% diag(sqrt(eigvals)))
# 
# # Check the results
# Z
# 
# # compare with pca
# pca_result <- prcomp(X)  # Standard PCA
# pca_scores <- pca_result$x
# 
# # Compare first few principal components
# print(head(Z)[1:6,1:4])  # Kernel PCA components
# print(head(pca_scores))  # Standard PCA components
# 
# # Plot the first two components for Standard PCA
# ggplot(pca_scores, aes(x = PC1, y = PC2, color = iris$Species)) +
#   geom_point(size = 3) +
#   labs(title = "Standard PCA (First Two Components)", x = "PC1", y = "PC2") +
#   scale_color_manual(values = c("setosa" = "red", "versicolor" = "green", "virginica" = "blue")) +
#   theme_minimal()
# 
# # Plot the first two components for Kernel PCA
# ggplot(Z, aes(x = V1, y = V2, color = iris$Species)) +
#   geom_point(size = 3) +
#   labs(title = "Kernel PCA (First Two Components)", x = "PC1", y = "PC2") +
#   scale_color_manual(values = c("setosa" = "red", "versicolor" = "green", "virginica" = "blue")) +
#   theme_minimal()

kpca <- function(data = data, response = response, kernel = c("linear", "gaussian", "diffusion"), numerical_sensitivity = 1e-10, gamma = NULL, tau = 1){
  
  # Set gamma based on kernel type if it's not provided
  if (is.null(gamma)) {
    if (kernel == "gaussian") {
      gamma <- 1 / (ncol(data) - 1)  # Default for gaussian kernel
    }
    else if (kernel == "diffusion") {
      gamma <- 1 / (ncol(data) - 1)  # Default for gaussian kernel
    }
  }
  # just take features and make into matrix
  X = as.matrix(data[,-which(colnames(data)==response)])
  
  # initialize kernel matrix
  # N = nrow(X)
  i = nrow(X) # individuals
  j = ncol(X) # genes
  K = matrix(data = 0, nrow = i, ncol = j)
  
  # Compute the RBF kernel: K(i, j) = exp(-gamma * ||xi - xj||^2)
  # for (i in 1:N) {
  #   for (j in 1:N) {
  #     dist_sq <- sum((X[i, ] - X[j, ])^2)
  #     K[i, j] <- exp(-gamma * dist_sq)
  #   }
  # }
  
  if(kernel == "linear"){
  # linear kernel for regular PCA
  K <- X %*% t(X)
  }
  else if(kernel == "gaussian"){
    K <- exp(-gamma * as.matrix(dist(X))^2)
  }
  else if (kernel == "diffusion") {
    
    L = 
    
    K = exp(beta*L)
  }
  
  
  # center kernel
  N <- nrow(K)
  H <- diag(N) - (1/N)*rep(1,N)%*%t(rep(1,N)) # centering matrix
  K_centered <- H%*%K%*%H # removes row means (HK) and col means (KH)
  
  # eigen decomposition
  eig <- eigen(K_centered)  # Compute eigenvalues and eigenvectors
  eigvals <- eig$values  # Eigenvalues
  eigvecs <- eig$vectors  # Eigenvectors
  
  # Set negative eigenvalues or very small values (close to zero) to zero
  eigvals[eigvals < numerical_sensitivity] <- 0  # Threshold for small eigenvalues (tune as needed)
  
  # Normalize eigenvectors
  Z <- as.data.frame(eigvecs %*% diag(sqrt(eigvals)))
  
  return(Z)
}

kpca_res <- kpca(data = iris, response = "Species", kernel = "diffusion", gamma = NULL, tau = 5)
# Plot the first two components for Kernel PCA
ggplot(kpca_res, aes(x = V1, y = V2, color = iris$Species)) +
  geom_point(size = 3) +
  labs(title = "Kernel PCA (First Two Components)", x = "PC1", y = "PC2") +
  scale_color_manual(values = c("setosa" = "red", "versicolor" = "green", "virginica" = "blue")) +
  theme_minimal()
