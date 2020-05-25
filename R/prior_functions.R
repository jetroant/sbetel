
#Prior generating function for var without constant
prior_fun_var <- function(model) {
  
  nn = floor(nrow(model$args$xy$yy)/model$args$p)
  xy <- model$args$xy
  p <- model$args$p
  m <- ncol(model$y)
  stat <- model$args$stat
  if(is.null(model$args$lambda)) lambda <- 1 else lambda <- model$args$lambda
  
  epsilon <- rep(NA, m)
  for(i in 1:m) epsilon[i] <- sqrt(ar(xy$yy[,i], aic = F, order.max = 1)$var.pred)
  
  prior_mean <- matrix(0, nrow = m*p, ncol = m)
  prior_mean[c(1:m),] <- diag(model$args$stat)
  
  xx <- xy$xx
  xx0 <- matrix(0, nrow = nn*p, ncol = m*p)
  for(i in 1:p) {
    rows <- (i*nn+1):(i*nn+nn)-nn
    cols <- (i*m+1):(i*m+m)-m
    rows_xx <- sample.int(nrow(xx), nn, replace = TRUE)
    xx0[rows,cols] <- xx[rows_xx, cols + 1]*i*lambda
  }
  
  yy0 <- xx0 %*% prior_mean
  for(i in 1:p) {
    rows <- (i*nn+1):(i*nn+nn)-nn
    cols <- (i*m+1):(i*m+m)-m
    for(j in 1:m) {
      yy0[rows,j] <- yy0[rows,j] + rnorm(nn, 0, epsilon[j])
    }
  }
  
  list("yy" = yy0, "xx" = xx0)
}

#Prior generating function for svar without constant
prior_fun_svar <- function(model) {
  
  if(is.null(model$args$nn)) {
    nn <- floor(nrow(model$args$xy$yy)/model$args$p)
  } else {
    nn <- model$args$nn
  }
  if(is.null(model$args$lambda)) {
    lambda <- 1
  } else {
    lambda <- model$args$lambda
  }
  if(is.null(model$args$prior_skew)) {
    prior_skew <- 0
  } else {
    prior_skew <- model$args$prior_skew
  }
  if(is.null(model$args$prior_dof)) {
    prior_dof <- 10
  } else {
    prior_dof <- model$args$prior_dof
  }
  xy <- model$args$xy
  p <- model$args$p
  m <- ncol(model$y)
  stat <- model$args$stat

  prior_mean <- matrix(0, nrow = m*p, ncol = m)
  prior_mean[c(1:m),] <- diag(model$args$stat)
  
  xx <- xy$xx
  xx0 <- matrix(0, nrow = nn*p, ncol = m*p)
  for(i in 1:p) {
    rows <- (i*nn+1):(i*nn+nn)-nn
    cols <- (i*m+1):(i*m+m)-m
    rows_xx <- sample.int(nrow(xx), nn, replace = TRUE)
    xx0[rows,cols] <- xx[rows_xx, cols + 1]*i*lambda
  }
  
  errors <- matrix(sgt::rsgt(nrow(xx0)*m, mu = 0, sigma = 1, 
                             lambda = prior_skew, p = 2, q = prior_dof/2,
                             mean.cent = TRUE, var.adj = TRUE),
                   ncol = m)
  errors <- errors %*% t(model$args$B0)
  yy0 <- xx0 %*% prior_mean + errors
  
  list("yy" = yy0, "xx" = xx0)
}

