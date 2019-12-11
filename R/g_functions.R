
g_var <- function(th, y, p, lambda = Inf, stat = NULL) {
  
  xy <- build_xy(y = y, p = p, lambda = lambda, stat = stat)
  n <- ncol(y)
  xx <- xy$xx
  yy <- xy$yy
  
  #AR block
  block1 <- matrix(NA, ncol = (n^2)*p+n, nrow = nrow(xx))
  u <- matrix(NA, ncol = n, nrow = nrow(xx))
  for(i in 1:n) {
    th_indices <- (i*(n*p+1)-(n*p+1)+1):(i*(n*p+1))
    u[,i] <- as.numeric(yy[,i] - xx %*% th[th_indices])
    block1[,th_indices] <- xx * u[,i]
  }
  
  #If the number of parameters equal the number of AR coefficients,
  #only moment conditions identifying those are returned.
  if(length(th) == ncol(y)*ncol(xx)) {
    gs <- block1
    return(gs)
  }
  
  #Sigma block
  B <- diag(n)
  B_indices <- (th_indices[length(th_indices)] + 1):length(th)
  B[!upper.tri(B)] <- th[B_indices]
  B_inv <- solve(B)
  errors <- u %*% t(B_inv)
  block2 <- matrix(NA, ncol = (n*(n+1)/2), nrow = nrow(errors))
  count <- 1
  for(i in 1:ncol(errors)) {
    for(j in 1:ncol(errors)) {
      if(i == j) {
        block2[,count] <- errors[,i] * errors[,j] - 1
        count <- count + 1
      }
      if(i > j) {
        block2[,count] <- errors[,i] * errors[,j]
        count <- count + 1
      }
    }
  }
  
  gs <- cbind(block1, block2)
  return(gs)
}
