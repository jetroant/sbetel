
initial_var <- function(xy) {
  xx <- xy$xx
  yy <- xy$yy
  ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  sigma <- t(yy - xx %*% ols_est) %*% (yy - xx %*% ols_est)/ nrow(yy)
  cross_xx_inv <- chol2inv( chol (crossprod(xx)))
  ols_cov <- kronecker(sigma, cross_xx_inv) 
  list("th" = c(ols_est), "cov" = ols_cov)
}
