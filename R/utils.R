
#Builds xx and yy from y and p (VAR)
build_xy_var <- function(y, p) {
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
  xx <- xx[-c(1:p),]
  yy <- y[-c(1:p),]
  ret <- list("xx" = xx, "yy" = yy)
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

#Load backup folder
load_backup <- function(dir) {
  
  #Load the backup directory
  wd <- getwd()
  setwd(dir)
  subchains <- list()
  filenames <- list.files()
  for(i in 1:length(filenames)) {
    subchains[[i]] <- readRDS(filenames[i])
  }
  setwd(wd)
  
  #Set parameters
  chain_length <- nrow(subchains[[1]]$chain)
  chain_number <- length(subchains)
  
  #Collect the chains and likelihoods
  allchains <- matrix(NA, 
                      ncol = ncol(subchains[[1]]$chain),
                      nrow = chain_length*chain_number
  )
  likelihoods <- c()
  for(i in 1:length(subchains)) {
    rows <- (i*chain_length+1):((i+1)*chain_length)-chain_length
    allchains[rows,] <- subchains[[i]]$chain
    burn <- length(subchains[[i]]$likelihoods) - nrow(subchains[[i]]$chain)
    if(burn > 0) {
      likelihoods <- c(likelihoods, subchains[[i]]$likelihoods[-c(1:burn)])
    } else {
      likelihoods <- c(likelihoods, subchains[[i]]$likelihoods)
    }
  }
  
  #Average acceptance rate
  accrates <- rep(NA, chain_number)
  for(i in 1:chain_number) accrates[i] <- subchains[[i]]$accrate
  
  #Collect the output
  ret <- list(sample = allchains,
              subchains = subchains,
              accrates = accrates,
              likelihoods = likelihoods,
              totaltime = NA,
              tune = NA,
              itermax = NA,
              type = NA)
  ret
}

#Approximates the prior (Make this parallel)
approx_prior <- function(model, nc = 100, cl = 1000, skip = FALSE) {
  
  st <- Sys.time()
  
  m <- ncol(model$y)
  p <- model$args$p
  prior_mean <- matrix(0, nrow = m*p, ncol = m)
  prior_mean[c(1:m),] <- diag(model$args$stat)
  prior_mean <- c(prior_mean, model$args$B0)
  
  model$args$constant <- FALSE
  og_nn <- model$args$nn
  
  if(skip == FALSE) {
    
    #Test run
    trys <- 3
    success <- FALSE
    while(trys > 0) {
      
      xy0 <- model$prior_fun(model)
      init <- tryCatch({
        sbetel:::initial_svar(xy0, model$args)
      }, error = function(e) {
        cat("GMM failed. Adjusting the sample size... \n")
        NULL
      })
      if(!is.null(init)) break
      trys <- trys - 1
      
      if(trys > 0) {
        model$args$nn <- model$args$nn*2
      }
      
    }
    if(is.null(init)) stop("GMM keeps failing. Prior could not be approximated.")
    if(!is.null(init)) cat("GMM succesfully estimated.")
  } 
  
  correction <- model$args$nn/og_nn
  
  cat("Sampling from asymptotic approximation of the prior... \n")
  pb <- txtProgressBar(min = 0, max = nc, style = 3)
  for(i in 1:nc) { 
    
    #Estimate new GMM model
    if(i > 1) {
      xy0 <- model$prior_fun(model)
      init <- tryCatch({
        sbetel:::initial_svar(xy0, model$args)
      }, error = function(e) {
        stop("GMM failed. \n")
      })
    }
    
    if(i == 1) {
      mat <- mvtnorm::rmvnorm(cl, 
                              mean = init$th, 
                              sigma = init$cov*correction)
    } else {
      temp <- mvtnorm::rmvnorm(cl, 
                               mean = init$th, 
                               sigma = init$cov*correction)
      mat <- rbind(mat, temp)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  time <- Sys.time() - st
  cat("Time took approximating the prior: ", round(time, 2), 
      " ", attributes(time)$units, "\n")
  mat
}









