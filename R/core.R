
#' @title Initialize a smoothed BETEL model
#' @description Initialize a smoothed (VAR) BETEL model by providing the data and 
#' the desired parameters. This function returns a list that defines
#' the model and can be passed either to \code{eval_sbetel()} or 
#' \code{est_sbetel()}.
#' 
#' @param y \code{(T x n)} numerical matrix with \code{n} time series
#' in its columns.
#' @param p Lag length. Defaults to \code{"auto"} in which case the optimal
#' lag length is chosen according to the out-of-sample mean squared 
#' forecasting error of the joint posterior mode. 
#' @param lambda The inverse amount of shrinkage towards the Minnesotaprior \code{(> 0)}. 
#' Defaults to \code{"auto"} in which case the optimal shrinkage is chosen according to the 
#' out-of-sample mean squared forecasting error of the joint posterior mode. 
#' No shrinkage is obtained by setting this to \code{Inf}.
#' @param bw Integer \code{(>= 0)}. The bandwidth parameter used for smoothing of the moment 
#' conditions. Defaults to \code{"auto"} in which case the optimal value is 
#' chosen according to \code{bwAndrews()} from \code{sandwich} package.
#' @param stat A numerical vector of length \code{n} \code{(0 <= n < 1)}. 
#' The prior mean of the coefficients for own first lags in the VAR model. 
#' For series with a priori more persistent dynamics, values closer to zero 
#' can be chosen. Defaults to zero for all series.
#' @param additional_priors \code{NULL} or numerical vector of length 2, where
#' both elements \code{> 0}. First element controls the inverse weight of
#' \emph{sum-of-coefficients} (SOC) prior whereas the other controls the 
#' \emph{initial-dummy-observation} (IDO) prior. Typical values in the literature 
#' for these values would be \code{c(1,1)}. Defaults to \code{NULL} in which case 
#' these priors are not used.
#' @return \code{init_sbetel()} returns a list that defines the model and can be passed
#' to \code{eval_sbetel()} or \code{est_sbetel}.
#' @examples
#' model <- init_sbetel(y)
#' @export
init_sbetel <- function(y,
                        p = "auto", 
                        lambda = "auto", 
                        bw = "auto",
                        stat = rep(0, ncol(y)), 
                        additional_priors = NULL) {
  
  #Chooses lags ('p') and shrinkage ('lambda')
  if(p == "auto") {
    p_grid <- 1:floor(nrow(y)/8)
    if(max(p_grid) > 13) p_grid <- 1:13
    if(min(p_grid) == 0) p_grid <- 1
  } else {
    p_grid <- p
  }
  fixed <- NULL
  if(lambda != "auto") fixed <- lambda
  objective <- Inf
  for(i in p_grid) {
    opt <- shrinkage_selector(y = y, 
                              p = i, 
                              stat = stat, 
                              additional_priors = additional_priors,
                              fixed = fixed)
    if(opt$objective < objective) {
      objective <- opt$objective
      lambda <- opt$minimum
      p <- i
    }
  }
  
  #OLS estimates as initial parameter values
  xy <- build_xy(y, p = p, lambda = lambda, stat = stat)
  xx <- xy$xx
  yy <- xy$yy
  OLS_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
  Sigma <- t(yy - xx %*% OLS_est) %*% (yy - xx %*% OLS_est)/ nrow(yy)
  cross_xx_inv <- chol2inv( chol (crossprod(xx)))
  OLS_cov <- kronecker(Sigma, cross_xx_inv)
  th_initial <- c(OLS_est, t(chol(Sigma))[!upper.tri(Sigma)])
  
  #Parameter covariance matrix from GMM for RWMH algorithm to use
  g_gmm <- function(th, x) {
    g_var(th = th, y = y, p = p, lambda = lambda, stat = stat)
  }
  gmm_model <- gmm::gmm(g_gmm, yy, t0 = th_initial, type = "twoStep", 
                        wmatrix = "ident", optfct = "nlminb")
  cov_initial <- gmm_model$vcov
  
  #Chooses the bandwidth for smoothing
  if(bw == "auto") {
    bw <- sandwich::bwAndrews(gmm_model, kernel = "Bartlett", prewhite = 0)
    bw <- floor(bw/2)
  }
  
  #Collects the model parameters etc.
  model <- list(y = y,
                g = g_var,
                p = p,
                lambda = lambda,
                bw = bw,
                stat = stat,
                additional_priors = additional_priors,
                xy = xy,
                th_initial = th_initial,
                cov_initial = cov_initial,
                type = "var")
  
  model
}

#' @title Evaluate a smoothed BETEL posterior density
#' @description Evaluate a (unnormalized) smoothed BETEL posterior density at parameter
#' values \code{th} given \code{model} (list) constructed by \code{init_sbetel()}.
#' 
#' @param th A numerical vector of parameter values.
#' @param model A list returned by \code{init_sbetel()}, defining the model.
#' @param itermax Maximum number of Newton-Rhapson iterations within the evaluation of the likelihood.
#' @return \code{eval_sbetel()} returns a numerical value that is the unnormalized
#' posterior density at \code{th}.
#' @examples
#' eval_sbetel(model$th_initial, model)
#' @export
eval_sbetel <- function(th, model, itermax = 20) {
  #Warnings are suppressed for now, 
  #as inv_sympd() from armadillo (c++) produces false warnings
  suppressWarnings(
    etel_rcpp(th = th, g = model$g, p = model$p, y = model$y, 
              bw = model$bw, lambda = model$lambda, td = model$xy$td,
              itermax = itermax)
  )
}

#Estimates the parameter posterior distribution using RWMH algorithm
#' @title Estimate a smoothed BETEL model
#' @description Estimates a smoothed BETEL model by generating a sample
#' from parameter posterior distribution by a RWMH algorithm.
#' 
#' @param model A list returned by \code{init_sbetel()}, defining the model.
#' @param N Number of draws from the posterior distribution.
#' @param tune A parameter that controls the step size of the RWMH algorithm.
#' Defaults to \code{"auto"} in which case sufficiently efficient parameter value
#' is automatically chosen by mimicking the manual tuning process with \code{N = 100}.
#' Causion is advised when using this default setting.
#' @param itermax Maximum number of Newton-Rhapson iterations within the evaluation 
#' of the likelihood.
#' @param burn \code{NULL} or numerical value \code{> 0}. The length of the 
#' burn-in sample. Defaults to \code{NULL} in which case there is no burn-in.
#' @param verbose Logical. If \code{TRUE} messages are produced. Defaults to \code{TRUE}. 
#' @return \code{est_sbetel()} returns a list containing the output of the 
#' RWMH algorithm.
#' @examples
#' model_output <- est_sbetel(model, N = 1000)
#' @export
est_sbetel <- function(model,
                       N,
                       tune = "auto",
                       itermax = 20, 
                       burn = NULL,
                       verbose = TRUE) {
  
  #Choosing the tuning parameter for suitable step size
  if(tune == "auto") {
    autotune_output <- auto_tune(model, verbose = verbose)
    tune <- autotune_output$tune
  }
  
  starttime <- Sys.time()
  
  #Adding burn-in
  if(!is.null(burn)) N <- N + burn
  
  #Pre-drawing the RW-steps
  if(verbose == TRUE) cat(paste0("Pre-drawing rw-steps... (N = ", N, ") \n"))
  moves <- mvtnorm::rmvnorm(N, mean = rep(0, length(model$th_initial)), sigma = model$cov_initial*tune)

  #Initializing RWMH chain
  mat <- matrix(NA, ncol = length(model$th_initial), nrow = N + 1)
  mat[1,] <- model$th_initial
  likelihoods <- rep(NA, nrow(mat))
  last_density <- eval_sbetel(model$th_initial, model, itermax)
  likelihoods[1] <- last_density
  
  #RWMH chain starts here
  if(verbose == TRUE) cat("Chain started, sampling from the posterior... \n")
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    proposal <- mat[i,] + moves[i,]
    proposal_density <- eval_sbetel(proposal, model, itermax)
    
    if(proposal_density - last_density > 0) {
      mat[i+1,] <- proposal
      last_density <- proposal_density
      
    } else if(proposal_density - last_density > log(runif(1, 0, 1))) {
      mat[i+1,] <- proposal
      last_density <- proposal_density
      
    } else {
      mat[i+1,] <- mat[i,]
    }
    
    if(!is.null(burn)) {
      if(i == burn) print("Burn-in complete...")
    }
    
    likelihoods[i+1] <- last_density
    if(verbose == TRUE) setTxtProgressBar(pb, i)
  }
  if(verbose == TRUE) close(pb)
  
  #Set colnames for the sample
  if(model$type == "var") {
    colnames(mat) <- c(paste0("A_", 1:(ncol(model$y)*ncol(model$xy$xx))),
                       paste0("C_", 1:((ncol(model$y)*(ncol(model$y)+1))/2))
    )
  }
  
  #Scrapping the burn-in sample
  if(!is.null(burn)) {
    mat_burned <- mat[c(1:burn),]
    likelihoods_burned <- likelihoods[c(1:burn)]
    mat <- mat[-c(1:burn),]
    likelihoods <- likelihoods[-c(1:burn)]
  } else {
    mat_burned <- NULL
    likelihoods_burned <- NULL
  }
  
  #Collecting everything
  accrate <- (length(unique(mat[,1]))-1)/(nrow(mat)-1)
  time <- Sys.time()-starttime
  output <- list(sample = NULL,
                 raw_sample = mat,
                 likelihoods = likelihoods,
                 accrate = accrate,
                 burned = list(sample = mat_burned, likelihoods = likelihoods_burned),
                 moves = moves,
                 time = time)
  if(tune == "auto") output$autotune <- autotune_output
  
  #Rescale the sample (if bw > 0)
  if(model$bw == 0) {
    output$sample <- output$raw_sample
    output$raw_sample <- NULL
  } else {
    beta_root <- sqrt(1/(model$bw*2+1))
    post_mean <- apply(mat, 2, mean)
    rescale <- function(row) (row - post_mean)*beta_root^(-1) + post_mean
    output$sample <- t(apply(mat, 1, rescale))
  }

  if(verbose == TRUE) cat(paste0("Acceptance rate: ", accrate, "\n"))
  if(verbose == TRUE) cat(paste0("Time: ", round(time, 3), " ", attributes(time)$units, "\n"))
  
  output
}

