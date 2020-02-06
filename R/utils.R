
#Builds X and Y from y and p (VAR)
build_xy <- function(y, p, lambda = Inf, stat = NULL, additional_priors = NULL) {
  
  #Build the xx and yy matrix
  n <- ncol(y)
  for(i in 1:p) {
    if(i == 1) {
      temp <- rbind(rep(NA, n), y[-(nrow(y)),])
      xx <- temp
    } else {
      temp <- rbind(rep(NA, n), temp[-(nrow(y)),])
      xx <- cbind(xx, temp)
    }
  }
  xx <- cbind(rep(1, nrow(xx)), xx)
  y0 <- y[c(1:p),]
  xx <- xx[-c(1:p),]
  yy <- y[-c(1:p),]
  
  #Minnesota prior via dummy observations
  if(lambda != Inf) {
    
    arsigmas <- c()
    for(i in 1:ncol(yy)) {
      armodel <- ar(yy[,i], aic = F, order.max = p)
      arsigmas <- c(arsigmas, sd(na.omit(armodel$resid)))
    }
    
    M0 <- rbind(diag(n), matrix(0, ncol = n, nrow = (n*p - n)))
    M0[which(M0 == 1)][stat] <- 0
    Psi0 <- diag(arsigmas)
    M1 <- lambda^(-1) * kronecker(diag(1:p), Psi0)
    
    td <- n*p+n+1
    yd <- rbind(lambda^(-1) * M0 %*% Psi0,
                sqrt(td)*Psi0,
                rep(0, n))
    xd <- rbind(cbind(rep(0, nrow(M1)), M1),
                matrix(0, ncol = ncol(M1) + 1, nrow = n),
                c(1e-6, rep(0, ncol(M1))))
    
    yy <- rbind(yd, yy)
    xx <- rbind(xd, xx)
  }
  
  #Additional priors
  if(!is.null(additional_priors)) {
    
    y_bar <- apply(y0, 2, mean)
    
    #SOC
    y_soc <- diag(y_bar/additional_priors[1])
    y_soc[stat,] <- 0
    x_soc <- cbind(0, y_soc)
    if(p > 1) for(i in 2:p) x_soc <- cbind(x_soc, y_soc)
    
    #IDO
    y_ido <- y_bar/additional_priors[2]
    x_ido <- c(1, rep(y_bar, p))/additional_priors[2]
    
    yy <- rbind(y_soc, y_ido, yy)
    xx <- rbind(x_soc, x_ido, xx)
    
    td <- n*p+n+1 + nrow(y_soc) + 1
  }
  
  ret <- list()
  ret$xx <- xx
  ret$yy <- yy
  if(lambda != Inf) ret$td <- td
  
  ret
}

#Returns all permutations of 'n' integers
permutations <- function(n){
  if(n == 1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow = n*p, ncol = n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i, sp + (sp >= i))
    }
    return(A)
  }
}

#Stacks VAR(p) coefficient matrix to VAR(1) coefficient matrix
stackA <- function(A) {
  A <- t(A)[,-1]
  m <- nrow(A)
  lags <- ncol(A)/m
  eye <- diag(m*lags-m)
  A <- rbind(A, cbind(eye, matrix(0, ncol = m, nrow= nrow(eye))))
  A
}

#For picking the time indices for narrative restrictions
pick_indices <- function(y_ts) cbind(y_ts, 1:nrow(y_ts))[,-c(1:ncol(y_ts))]



