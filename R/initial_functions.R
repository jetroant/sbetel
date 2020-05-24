
initial_var <- function(xy, args = NULL) {
  xx <- xy$xx
  yy <- xy$yy
  ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  sigma <- t(yy - xx %*% ols_est) %*% (yy - xx %*% ols_est)/ nrow(yy)
  cross_xx_inv <- chol2inv( chol (crossprod(xx)))
  ols_cov <- kronecker(sigma, cross_xx_inv) 
  list("th" = c(ols_est), "cov" = ols_cov)
}

initial_svar <- function(xy, args, bw = FALSE) {
  
  #First stage initial values for the autoregressive parameters from OLS
  xx <- xy$xx
  yy <- xy$yy
  ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  sigma <- t(yy - xx %*% ols_est) %*% (yy - xx %*% ols_est)/ nrow(yy)
  cross_xx_inv <- chol2inv( chol (crossprod(xx)))
  ols_cov <- kronecker(sigma, cross_xx_inv) 
  ols <- list("th" = c(ols_est), "cov" = ols_cov)
  init_B <- diag(args$epsilon)
  init_th <- c(ols$th, init_B)
  
  #Second stage initial values from GMM
  args$xy <- xy
  g_gmm <- function(th, x) {
    sbetel:::g_svar(th = th, y = x, args = args)
  }
  gmm_init <- gmm::gmm(g_gmm, yy, t0 = init_th, type = "twoStep", 
                       wmatrix = "optimal", optfct = "nlminb")
  ret <- list("th" = gmm_init$coefficients, "cov" = gmm_init$vcov)
  
  if(bw == TRUE) {
    bw <- sandwich::bwAndrews(gmm_init, kernel = "Bartlett", prewhite = 0)
    bw <- floor(bw/2)
    ret$bw <- bw
  }
  
  ret
}


