
g_var <- function(th, y, args) {
  
  p <- args$p
  lambda <- args$lambda
  stat <- args$stat
  
  xy <- build_xy(y = y, p = p, lambda = lambda, stat = stat)
  n <- ncol(y)
  xx <- xy$xx
  yy <- xy$yy
  
  block2 <- NULL
  block3 <- NULL
  block4 <- NULL
  
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
  if(length(th) > ncol(y)*ncol(xx)) { 
    
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
    
  } 
  
  if(!is.null(args$moment_conditions)) {
    
    cMat <- args$moment_conditions
    
    if(2 %in% cMat) {
      prod2 <- function(row) prod(row[which(row != 0)])
      cMat2 <- cMat[which(apply(cMat, 1, prod2) == 4),]
      cMat <- cMat[which(apply(cMat, 1, prod2) == 3),]
      
      if(is.null(nrow(cMat2))) cMat2 <- matrix(cMat2, nrow = 1)
      if(is.null(nrow(cMat))) cMat <- matrix(cMat, nrow = 1)
      
      block4 <- matrix(NA, ncol = nrow(cMat2), nrow = nrow(errors))
      for(i in 1:nrow(cMat2)) {
        for(j in 1:ncol(cMat2)) {
          if(j == 1) {
            cond <- errors[,j]^cMat2[i,j]
          } else {
            cond <- cond * errors[,j]^cMat2[i,j]
          }
        }
        block4[,i] <- cond - 1
      }
    }
    
    block3 <- matrix(NA, ncol = nrow(cMat), nrow = nrow(errors))
    for(i in 1:nrow(cMat)) {
      for(j in 1:ncol(cMat)) {
        if(j == 1) {
          cond <- errors[,j]^cMat[i,j]
        } else {
          cond <- cond * errors[,j]^cMat[i,j]
        }
      }
      block3[,i] <- cond
    }
    
  }
  
  gs <- cbind(block1, block2, block3, block4)
  return(gs)
}
