
g_var <- function(th, y, args) {
  
  n <- ncol(y)
  p <- args$p
  xx <- args$xy$xx
  yy <- args$xy$yy
  
  if(args$constant == TRUE) {
    block1 <- matrix(NA, ncol = (n^2)*p+n, nrow = nrow(xx))
    u <- matrix(NA, ncol = n, nrow = nrow(xx))
    for(i in 1:n) {
      th_indices <- (i*(n*p+1)-(n*p+1)+1):(i*(n*p+1))
      u[,i] <- as.numeric(yy[,i] - xx %*% th[th_indices])
      block1[,th_indices] <- xx * u[,i]
    }
    
  } else {
    block1 <- matrix(NA, ncol = (n^2)*p, nrow = nrow(xx))
    u <- matrix(NA, ncol = n, nrow = nrow(xx))
    for(i in 1:n) {
      th_indices <- ((i*(n*p))-n*p+1):(i*(n*p))
      u[,i] <- as.numeric(yy[,i] - xx %*% th[th_indices])
      block1[,th_indices] <- xx * u[,i]
    }
  }

  block1
}

g_svar <- function(th, y, args) {
  
  n <- ncol(y)
  p <- args$p
  xx <- args$xy$xx
  yy <- args$xy$yy
  
  block3 <- NULL
  block4 <- NULL
  block5 <- NULL
  
  if(args$constant == TRUE) {
    block1 <- matrix(NA, ncol = (n^2)*p+n, nrow = nrow(xx))
    u <- matrix(NA, ncol = n, nrow = nrow(xx))
    for(i in 1:n) {
      th_indices <- (i*(n*p+1)-(n*p+1)+1):(i*(n*p+1))
      u[,i] <- as.numeric(yy[,i] - xx %*% th[th_indices])
      block1[,th_indices] <- xx * u[,i]
    }
    
  } else {
    block1 <- matrix(NA, ncol = (n^2)*p, nrow = nrow(xx))
    u <- matrix(NA, ncol = n, nrow = nrow(xx))
    for(i in 1:n) {
      th_indices <- ((i*(n*p))-n*p+1):(i*(n*p))
      u[,i] <- as.numeric(yy[,i] - xx %*% th[th_indices])
      block1[,th_indices] <- xx * u[,i]
    }
  }
  
  #B-matrix and errors
  if(length(th) <= th_indices[length(th_indices)]) stop("Initial parameter values misspecified.")
  B_indices <- (th_indices[length(th_indices)] + 1):length(th)
  B <- matrix(th[B_indices], ncol = n)
  B_inv <- solve(B)
  errors <- u %*% t(B_inv) 
  
  #Second order conditions
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
  
  #Third and fourth order moment conditions
  cMat <- sbetel:::permutations(n)
  if("skewness" %in% args$identification) {
    block3 <- matrix(NA, nrow = nrow(errors), ncol = nrow(cMat))
    for(i in 1:nrow(cMat)) {
      block3[,i] <- errors[,cMat[i,1]]^2 * errors[,cMat[i,2]]
    }
  }
  if("kurtosis" %in% args$identification) {
    block4 <- matrix(NA, nrow = nrow(errors), ncol = nrow(cMat))
    for(i in 1:nrow(cMat)) {
      block4[,i] <- errors[,cMat[i,1]]^3 * errors[,cMat[i,2]]
    }
  }
  if("symmetric" %in% args$identification) {
    cMat_sym <- matrix(NA, nrow = (n-1)*n/2, ncol = 2)
    count <- 0
    for(i in 1:n) {
      for(j in 1:n) {
        if(i > j) {
          count <- count + 1
          cMat_sym[count,] <- c(i,j)
        }
      }
    }
    block5 <- matrix(NA, nrow = nrow(errors), ncol = nrow(cMat_sym))
    for(i in 1:nrow(cMat_sym)) {
      block5[,i] <- errors[,cMat_sym[i,1]]^2 * errors[,cMat[i,2]]^2 - 1
    }
  }
  
  gs <- cbind(block1, block2, block3, block4, block5)
  return(gs)
}
