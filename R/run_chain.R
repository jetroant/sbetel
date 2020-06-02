
run_chain <- function(chain_name,
                      type,
                      model,
                      chain_length,
                      tune,
                      burn,
                      backup,
                      itermax,
                      wd = getwd(),
                      trys = 1,
                      progressbar = FALSE) {
  
  #Add burn-in
  N <- chain_length + burn 
  
  #If initial values from GMM fall outside the support of posterior
  #density, the prior data and initial values are regenerated
  success <- FALSE
  while(trys > 0) {
    trys <- trys - 1
    
    #Select type
    if(type == "posterior") {
      xy0 <- model$prior_fun(model)
      xy <- model$args$xy
      if(ncol(xy0$xx) != ncol(xy$xx)) {
        xy0$xx <- cbind(0, xy0$xx)
      }
      xy$xx <- rbind(xy0$xx, xy$xx)
      xy$yy <- rbind(xy0$yy, xy$yy)
      td <- nrow(xy0$yy)
      
    } else if(type == "prior") {
      xy <- model$prior_fun(model)
      model$bw <- 0
      model$args$constant <- FALSE
      td <- 0
      
    } else if(type == "likelihood") {
      xy <- model$args$xy
      td <- 0
    }
    
    #Initial values and scale
    init <- tryCatch({
      model$initial(xy, model$args)
    }, error = function(e) {
      mes <- paste0("Initial parameter values cannot be calculated \n",
                    "  (most likely GMM fails): \n",
                    e)
      stop(mes)
    })
    init <- model$initial(xy, model$args)
    initial_th <- init$th
    initial_scale <- init$cov
    
    last_density <- eval_sbetel(th = initial_th, 
                                model = model, 
                                xy = xy,
                                td = td,
                                itermax = itermax)
    
    if(last_density != -Inf) {
      success <- TRUE
      break
    }
  }
  if(success == FALSE) {
    stop("Initial parameter values keep falling outside the support of the posterior density.")
  }
  
  #Pre-draw RW-steps
  moves <- mvtnorm::rmvnorm(N, 
                            mean = rep(0, length(initial_th)), 
                            sigma = initial_scale)*tune
  
  #Initialize the chain
  mat <- matrix(NA, ncol = length(initial_th), nrow = N + 1)
  mat[1,] <- initial_th
  likelihoods <- rep(NA, nrow(mat))
  likelihoods[1] <- last_density
  
  #Start the chain
  if(progressbar == TRUE) pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    
    proposal <- mat[i,] + moves[i,]
    proposal_density <- eval_sbetel(th = proposal, 
                                    model = model, 
                                    xy = xy,
                                    td = td,
                                    itermax = itermax)
    
    if(proposal_density - last_density > 0) {
      mat[i+1,] <- proposal
      last_density <- proposal_density
      
    } else if(proposal_density - last_density > log(runif(1, 0, 1))) {
      mat[i+1,] <- proposal
      last_density <- proposal_density
      
    } else {
      mat[i+1,] <- mat[i,]
    }
    
    likelihoods[i+1] <- last_density
    if(progressbar == TRUE) setTxtProgressBar(pb, i)
  }
  if(progressbar == TRUE) close(pb)
  
  #Rescale the sample (if bw > 0)
  if(model$bw != 0) {
    beta_root <- sqrt(1/(model$bw*2+1))
    post_mean <- apply(mat, 2, mean)
    rescale <- function(row) (row - post_mean)*beta_root^(-1) + post_mean
    mat <- t(apply(mat, 1, rescale))
  }
  
  #Burn-in
  burned <- mat[1:(burn+1),]
  mat <- mat[-c(1:(burn+1)),]
  
  #Acceptance rate
  accrate <- (length(unique(mat[,1]))-1)/(nrow(mat)-1)
  
  #Collect the output
  if(type != "posterior") xy0 <- NULL
  ret <- list(chain = mat,
              burned = burned,
              likelihoods = likelihoods,
              accrate = accrate,
              xy0 = xy0)
  
  #Save the output
  if(!is.null(backup)) {
    setwd(wd)
    saveRDS(ret, paste0(backup, "/chain_", chain_name, ".rds"))
  }
  
  ret
}