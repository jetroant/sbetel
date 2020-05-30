
initial_var <- function(xy, args = NULL) {
  xx <- xy$xx
  yy <- xy$yy
  ols_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  sigma <- t(yy - xx %*% ols_est) %*% (yy - xx %*% ols_est)/ nrow(yy)
  cross_xx_inv <- chol2inv( chol (crossprod(xx)))
  ols_cov <- kronecker(sigma, cross_xx_inv) 
  list("th" = c(ols_est), "cov" = ols_cov)
}

initial_svar <- function(xy, args, bw = FALSE, verbose = FALSE) {
  
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
  if(verbose == TRUE) trace <- 1 else trace <- 0
  if(args$wmatrix == "optimal") {
    gmm_init <- gmm:::gmm(g_gmm, yy, t0 = init_th, type = "twoStep", 
                         wmatrix = args$wmatrix, optfct = "nlminb",
                         control = list(trace = trace))
    
  } else {
    gmm_init <- gmm::gmm(g_gmm, yy, t0 = init_th, type = "twoStep", 
                         wmatrix = args$wmatrix, optfct = "nlminb",
                         prewhite = 0, kernel = "Bartlett",
                         control = list(trace = trace))
  }
  
  ret <- list("th" = gmm_init$coefficients, "cov" = gmm_init$vcov)
  if(bw == TRUE) {
    bw <- sandwich::bwAndrews(gmm_init, kernel = "Bartlett", prewhite = 0)
    ret$bw_raw <- bw
    ret$bw <- floor(bw/2)
  }
  ret
}

initial_svar_robust <- function(xy, args) {
  
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
  
  #Marginal optimization of B
  init_B <- diag(args$epsilon)
  diag_ind <- which(init_B != 0) 
  lower <- rep(-Inf, length(init_B))
  lower[diag_ind] <- 0.01
  obj <- function(x) {
    gs <- args$g(c(ols$th, x), yy, args)
    mvec <- apply(gs, 2, mean)
    t(mvec) %*% diag(length(mvec)) %*% mvec
  }
  opt <- nlminb(init_B, obj, lower = lower)
  opt_B <- matrix(opt$par, ncol = ncol(yy))
  init_th <- c(ols$th, opt_B)
  
  #Approximative block-covariance matrix
  gs_B <- args$g(init_th, yy, args)[,-c(1:ncol(ols_cov))] / nrow(yy)
  grad <- function(x) {
    gs_B <- args$g(c(ols$th, x), yy, args)[,-c(1:ncol(ols_cov))]
    mvec <- apply(gs_B, 2, mean)
  }
  G <- numDeriv::jacobian(grad, opt_B)
  B_cov <- chol2inv(chol( t(G) %*% chol2inv(chol(crossprod(gs_B))) %*% G ))
  init_cov <- rbind(cbind(ols_cov, matrix(0, ncol = ncol(B_cov), nrow = nrow(ols_cov))),
                    cbind(matrix(0, ncol = ncol(ols_cov), nrow = nrow(B_cov)), B_cov))
  
  ret <- list("th" = init_th, "cov" = init_cov)
  ret
}



