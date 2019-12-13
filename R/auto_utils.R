
#Selects optimal shrinkage for BVAR with minnesota prior
shrinkage_selector <- function(y, p, stat = rep(0, ncol(y)), additional_priors = NULL, fixed = NULL) {
  
  toMin <- function(lambda) {
    xy <- build_xy(y, p, lambda = lambda, stat = stat, additional_priors = additional_priors)
    td <- ifelse(is.null(xy$td), 0, xy$td)
    xx <- xy$xx
    yy <- xy$yy
    errors <- matrix(NA, ncol = ncol(yy), nrow = nrow(yy))
    for(i in (1+xy$td):nrow(yy)) {
      OLS_est <- chol2inv(chol(crossprod(xx[-i,]))) %*% t(xx[-i,]) %*% yy[-i,]
      errors[i,] <-  xx[i,] %*% OLS_est - yy[i,]
    }
    if(xy$td != 0) errors <- errors[-c(1:xy$td),]; yy <- yy[-c(1:xy$td),]
    mse <- apply(errors^2, 2, mean)
    weights <- diag(apply(yy, 2, var)^(-1))
    mse %*% weights %*% mse
  }
  
  if(is.null(fixed)) {
    opt <- optimize(toMin, c(0.01, 10))
    return(opt)
  } else {
    return(list(minimum = fixed, objective = toMin(fixed)))
  }
}

#Approximates an efficient step size for RWMH algortihm,
#mimicing the manual tuning process
auto_tune <- function(model, target = c(0.22, 0.28), 
                      limit = 20, n = 100, verbose = TRUE) {
  if(verbose == TRUE) cat(paste0("Auto tuning step size... \n"))
  found <- FALSE
  tune <- 2.38^2/length(model$initial$th)
  step <- 0.1
  direction <- c()
  while(found == FALSE) {
    
    #Autotune fails if limit is reached
    if(length(direction) > limit) stop("Autotune failed. Limit reached.")
    
    s <- est_sbetel(model, N = n, tune = tune, verbose = FALSE)
    if(s$output$accrate > target[2]) {
      
      #Step size halves if the direction of the search changes
      direction <- c(1, direction)
      if(length(direction) != 1 & direction[1] == -direction[2]) step <- step/2
      tune <- tune + step
    } else if(s$output$accrate < target[1]) {
      
      #Step size halves if the direction of the search changes
      direction <- c(-1, direction)
      if(length(direction) != 1 & direction[1] == -direction[2]) step <- step/2
      if((tune - step) < 0) {
        tune <- tune/2
      } else {
        tune <- tune - step
      }
      
    } else {
      
      #Parameter is accepted if 'accrate' is still within the target
      s <- est_sbetel(model, N = 100, tune = tune, verbose = FALSE)
      if(s$output$accrate > target[1] & s$output$accrate < target[2]) {
        found <- TRUE
      }
    }
  }
  
  if(verbose == TRUE) cat(paste0("Autotune: 'tune' = ", tune, " / 'accrate' ~ ", round(s$output$accrate, 2), " / 'n' = ", n, "\n"))
  list(accrate = s$output$accrate, tune = tune, iters = length(direction))
}

