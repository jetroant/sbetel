
#Prior generating function for var without constant
prior_fun_var <- function(model, nn = 50, epsilon = 0.1, sigma_x = 0.1) {
  
  p <- model$args$p
  m <- ncol(model$y)
  stat <- model$args$stat
  
  prior_mean <- matrix(0, nrow = m*p, ncol = m)
  prior_mean[c(1:m),] <- diag(model$args$stat)
  
  xx0 <- matrix(0, nrow = nn*p, ncol = m*p)
  for(i in 1:p) {
    rows <- (i*nn+1):(i*nn+nn)-nn
    cols <- (i*m+1):(i*m+m)-m
    xx0[rows,cols] <- rnorm(nn*m, 0, sigma_x*i)
  }
  
  yy0 <- xx0 %*% prior_mean
  for(i in 1:p) {
    rows <- (i*nn+1):(i*nn+nn)-nn
    cols <- (i*m+1):(i*m+m)-m
    yy0[rows,] <- yy0[rows,] + rnorm(nn*m, 0, epsilon)
  }
  
  list("yy" = yy0, "xx" = xx0)
}

