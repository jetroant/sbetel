
#Estimate Gaussian-inverse-Wishart BVAR
est_gaussian_bvar <- function(model, N) {
  
  xx <- model$args$xy$xx
  yy <- model$args$xy$yy
  ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  u <- yy - xx %*% ols_est

  #Sample from inverse-Wishart
  scale <- t(u) %*% (u)
  sample_P <- matrix(NA, ncol = length(scale), nrow = N)
  for(i in 1:N) sample_P[i,] <- chol(MCMCpack::riwish(v = nrow(u), scale))
  
  #Sample from conditional multinormal
  sample_A <- matrix(NA, ncol = length(ols_est), nrow = N)
  cross_xx_inv <- chol2inv(chol(crossprod(xx)))
  for(i in 1:N) {
    sampled_P <- matrix(sample_P[i,], ncol = ncol(yy))
    sampled_sigma <- t(sampled_P) %*% sampled_P
    x_vec <- rnorm(n = length(ols_est), 0, 1)
    sample_A[i,] <- c(ols_est) + t(kronecker(chol(sampled_sigma), chol(cross_xx_inv))) %*% x_vec
  }
  
  sample_P <- sample_P[,which(sample_P[1,] != 0)]
  gaussian_sample <- cbind(sample_A, sample_P)
  colnames(gaussian_sample) <- c(paste0("A_", 1:ncol(sample_A)),
                                 paste0("P_", 1:ncol(sample_P)))
  
  model$output$gaussian_sample <- gaussian_sample
  model
}
