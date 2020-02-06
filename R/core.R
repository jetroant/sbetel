#'@useDynLib sbetel, .registration = TRUE
#'@importFrom Rcpp evalCpp

#' @title Initialize a smoothed BETEL model
#' @description Initialize a smoothed BETEL model by providing the model defining 
#' function \code{g()} and the data \code{y} to be used for estimation of the model. Alternatively 
#' ready-made functions can be used (e.g. by setting \code{g = "var"}). \code{init_sbetel()} returns 
#' a list that defines the model and that can be passed on either to \code{eval_sbetel()} or 
#' \code{est_sbetel()}.
#' 
#' @param g A function of the form \code{g(th, y, args)} defining the model.
#' Inputs should include a numerical parameter vector \code{th}, a data matrix 
#' \code{y} and optional further arguments through a list called \code{args}.
#' Output shoud be a numerical matrix with its columns corresponding to the sample
#' moment conditions of the model. Defaults to \code{"var"} in which case a 
#' ready-made function for vector autoregressive models is used. For more, 
#' see details.
#' @param y A data matrix. Input for \code{g()} above. E.g. for \code{g = "var"} this
#' would be a \code{T x k} matrix with \code{T} observations of \code{k} variables.
#' @param bw Integer \code{(>= 0)}. The bandwidth parameter controlling smoothing of the moment 
#' conditions. Defaults to "auto" in which case the optimal parameter value is 
#' chosen according to \code{bwAndrews()}. For more, see details.
#' @param args A list containing optional further arguments to \code{g()}.
#' If \code{g = "var"} is chosen, this can be left undefined in which case the 
#' parameters of the var model are automatically chosen to minimize the 
#' out-of-sample squared forecasting errors of the model. For example, the lag 
#' length and the shrinkage parameter can be manually fixed by defining 
#' \code{args = list(p = 5, lambda = 0.2)}. For more, see details.
#' @param initial A list of initial parameter values that defaults to 
#' \code{list(th = NULL, cov = NULL)}. At least \code{initial$th} must be predefined, 
#' if \code{g} is not set as \code{"var"}. If \code{initial$cov} is not predefined, 
#' it is approximated with a two-step GMM.
#' @param verbose Logical. Defaults to \code{TRUE} in which case messages are produced.
#' @return \code{init_sbetel()} returns a list that defines a sbetel model 
#' and can be passed on to \code{eval_sbetel()} or \code{est_sbetel()}.
#' @details TBA
#' @examples
#'\dontrun{
#' model <- init_sbetel(g = "var", y = y)
#'}
#' @export
init_sbetel <- function(g = "var",
                        y,
                        bw = "auto",
                        args = list(sigma = TRUE,
                                    gmm = FALSE),
                        initial = list(th = NULL, 
                                       cov = NULL),
                        verbose = TRUE
                        ) {
  
  if(is.null(args$gmm)) args$gmm <- FALSE
  
  if(is.character(g)) {
    
    if(g != "var") stop("Only 'var' defined for 'g' of type 'character' in 'init_sbetel'")
    
    g <- g_var
    type <- "var"
    
    if(is.null(args$stat)) args$stat <- rep(0, ncol(y))
    if(is.null(args$sigma)) args$sigma <- FALSE
    if(is.null(args$additional_priors)) args$additional_priors <- NULL
    
    #Chooses lags ('p') and shrinkage ('lambda')
    if(is.null(args$p)) {
      p_grid <- 1:floor(nrow(y)/8)
      if(max(p_grid) > 13) p_grid <- 1:13
      if(min(p_grid) == 0) p_grid <- 1
    } else {
      p_grid <- p
    }
    fixed <- NULL
    if(!is.null(args$lambda)) fixed <- lambda
    objective <- Inf
    for(i in p_grid) {
      opt <- shrinkage_selector(y = y, 
                                p = i, 
                                stat = args$stat, 
                                additional_priors = args$additional_priors,
                                fixed = fixed)
      if(opt$objective < objective) {
        objective <- opt$objective
        lambda <- opt$minimum
        p <- i
      }
    }
    args$p <- p
    args$lambda <- lambda
    
    #OLS estimates as initial parameter values
    xy <- build_xy(y, p = args$p, lambda = args$lambda, stat = args$stat)
    args$xy <- xy
    xx <- xy$xx
    yy <- xy$yy
    OLS_est <- chol2inv(chol(crossprod(xx))) %*% t(xx) %*% yy
    Sigma <- t(yy - xx %*% OLS_est) %*% (yy - xx %*% OLS_est)/ nrow(yy)
    cross_xx_inv <- chol2inv( chol (crossprod(xx)))
    OLS_cov <- kronecker(Sigma, cross_xx_inv)
    if(args$sigma == TRUE & is.null(initial$th)) {
      initial$th <- c(OLS_est, t(chol(Sigma))[!upper.tri(Sigma)])
    } else if(is.null(initial$th)) {
      initial$th <- c(OLS_est)
    }
    
  } else {
    
    type <- "custom"
    if(is.null(initial$th)) {
      stop("Initial parameter values must be provided if g = 'var' is not used!")
    }
    
  }
  
  #Parameter covariance matrix from GMM or OLS+bootstrap for RWMH algorithm to use
  #(if not provided by user)
  if(is.null(initial$cov) | !is.null(args$moment_conditions)) {
    
    if(type != "var" | args$gmm == TRUE) {
      g_gmm <- function(th, x) {
        g(th = th, y = y, args = args)
      }
      if(verbose == TRUE) cat("Estimating GMM model for initial values... \n")
      gmm_model <- gmm::gmm(g_gmm, y, t0 = initial$th, type = "twoStep", 
                            wmatrix = "optimal", optfct = "nlminb", vcov = "HAC")
      initial$cov <- gmm_model$vcov
      if(args$gmm == TRUE) initial$th <- gmm_model$coefficients
      
    } else {
      boot_n <- 1000
      temp <- matrix(NA, nrow = boot_n, ncol = length(initial$th) - length(OLS_est))
      u <- yy - xx %*% OLS_est
      for(i in 1:boot_n) {
        u_booted <- u[sample.int(nrow(u), nrow(u), replace = TRUE),]
        Sigma_booted <- crossprod(u_booted) / nrow(u_booted)
        temp[i,] <- t(chol(Sigma_booted))[!upper.tri(Sigma_booted)]
      }
      P_cov <- diag(apply(temp, 2, var))
      initial$cov <- rbind(cbind(OLS_cov, matrix(0, ncol = ncol(P_cov), nrow = nrow(OLS_cov))),
                           cbind(matrix(0, ncol = ncol(OLS_cov), nrow = nrow(P_cov)), P_cov))
    }
  }
  
  #Check for computational singularity of the initial covariance matrix
  #and load the diagonal if necessary to obtain a non-singular initial 
  #covariance matrix estimator.
  if(min(eigen(initial$cov)$values) < 0) {
    while(min(eigen(initial$cov)$values) < 0) {
      diag(initial$cov) <- diag(initial$cov)*1.01
    }
    if(verbose == TRUE) cat("Diagonal of the initial parameter covariance matrix loaded to avoid computational singularity. \n")
  }
  
  #Chooses the bandwidth parameter for smoothing
  if(bw == "auto") {
    g_gmm <- function(th, x) {
      g(th = th, y = y, args = args)
    }
    if(verbose == TRUE) cat("Choosing the optimal bandwidth parameter... \n")
    gmm_fs <- gmm::gmm(g_gmm, y, t0 = initial$th, type = "twoStep", 
                       wmatrix = "ident", optfct = "nlminb")
    bw <- sandwich::bwAndrews(gmm_fs, kernel = "Bartlett", prewhite = 0)
    bw <- floor(bw/2)
  }
  
  #Collects the model parameters etc.
  model <- list(g = g,
                y = y,
                bw = bw,
                args = args,
                initial = initial,
                type = type)
  
  #Finally, especially in case of additional moment conditions, 
  #feasibility/optimality of the initial parameter values must be checked
  if(!is.null(args$moment_conditions)) {
    
    #Evaluate the likelihood at initial values
    init_like <- eval_sbetel(model$initial$th, model)
    
    #Potentially optimal initial values are searched for with a 'marginal gmm'
    initial_fixed <- initial$th[1:length(OLS_est)]
    initial_not_fixed <- initial$th[-c(1:length(OLS_est))]
    g_obj <- function(th) {
      gs <- g(th = c(initial_fixed, th), y = y, args = args)[,-c(1:length(initial$th))]
      gs_mean <- apply(gs, 2, mean)
      gs_mean %*% diag(length(gs_mean)) %*% gs_mean
    }
    opt <- nlminb(start = initial_not_fixed, objective = g_obj)
    new_initial <- c(initial_fixed, opt$par)
    
    #Evaluate new initial likelihood
    init_like_new <- eval_sbetel(new_initial, model)
    
    #Check the optimality of the initial values
    if(init_like_new > init_like) {
      model$initial$th <- new_initial
      if(verbose == TRUE) cat("Optimal initial parameter values found with 'marginal gmm'. \n")
    }
  }
  
  #Evaluate the likelihood at final initial values
  init_like <- eval_sbetel(model$initial$th, model)
  
  #Check the feasibility of the initial values
  if(init_like == -Inf) {
    stop("No feasible initial parameter values found. This suggests either that (i) the specified model is not supported by the data, (ii) there is not enough data for identification or (iii) the initial parameter values are poorly chosen.")
  }
  
  model
}

#' @title Evaluate a smoothed BETEL posterior density
#' @description Evaluate a (unnormalized) smoothed BETEL posterior density at parameter
#' values \code{th} given the \code{model} (list) generated by \code{init_sbetel()}. 
#' A wrapper for c++ implementation.
#' 
#' @param th A numerical vector of parameter values.
#' @param model A list returned by \code{init_sbetel()} defining the model.
#' @param itermax Maximum number of Newton-Rhapson iterations within the evaluation of the likelihood. 
#' Defaults to \code{itermax = 20}.
#' @return \code{eval_sbetel()} returns a numerical value that is the unnormalized
#' posterior density at \code{th}.
#' @examples
#'\dontrun{
#' eval_sbetel(th = model$initial$th, model = model)
#'}
#' @export
eval_sbetel <- function(th, model, itermax = 20) {
  if(model$type == "var") {
    td <- model$args$xy$td
  } else {
    td <- 0
  } 
  etel_rcpp(th = th, 
            g = model$g, 
            y = model$y, 
            bw = model$bw, 
            td = td,
            itermax = itermax,
            args = model$args
  )
}

#' @title Estimate a smoothed BETEL model
#' @description Estimate a smoothed BETEL model by generating a sample
#' from the parameter posterior distribution with a RWMH algorithm.
#' 
#' @param model A list returned by \code{init_sbetel()} defining the model.
#' @param N Number of draws from the posterior distribution.
#' @param tune A parameter that controls the step size of the RWMH algorithm.
#' Defaults to \code{"auto"} in which case sufficiently efficient parameter value
#' is automatically chosen by mimicking the manual tuning process with \code{N = n}.
#' Causion is advised when using this default setting.
#' @param n Only used if \code{tune = "auto"}. The number of draws used for auto tuning the 
#' step size. Defaults to \code{n = 100}.
#' @param itermax Maximum number of Newton-Rhapson iterations within the evaluation 
#' of the likelihood. Defaults to \code{itermax = 20}.
#' @param burn \code{NULL} or numerical value \code{> 0}. The length of the 
#' burn-in sample. Defaults to \code{NULL} in which case there is no burn-in.
#' @param verbose Logical. Defaults to \code{TRUE} in which case messages are produced.
#' @return \code{est_sbetel()} returns a list containing the model and the output of the 
#' RWMH algorithm.
#' @details TBA
#' @examples
#' \dontrun{
#' model_output <- est_sbetel(model = model, N = 1000)
#' }
#' @export
est_sbetel <- function(model,
                       N,
                       tune = "auto",
                       itermax = 20, 
                       burn = NULL,
                       n = 100,
                       verbose = TRUE,
                       backup = NULL,
                       chain_name = 0) {
  
  #Choosing the tuning parameter for suitable step size
  if(tune == "auto") {
    autotune_output <- auto_tune(model, n = n, verbose = verbose)
    tune <- autotune_output$tune
  }
  
  starttime <- Sys.time()
  
  #Adding burn-in
  if(!is.null(burn)) N <- N + burn
  
  #Pre-drawing the RW-steps
  if(verbose == TRUE) cat(paste0("Pre-drawing rw-steps... (N = ", N, ") \n"))
  moves <- mvtnorm::rmvnorm(N, mean = rep(0, length(model$initial$th)), sigma = model$initial$cov*tune)

  #Initializing RWMH chain
  mat <- matrix(NA, ncol = length(model$initial$th), nrow = N + 1)
  mat[1,] <- model$initial$th
  likelihoods <- rep(NA, nrow(mat))
  last_density <- eval_sbetel(model$initial$th, model, itermax)
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
    
    if(!is.null(backup)) {
      if(i %% backup == 0) saveRDS(list(mat, likelihoods), paste0("sbetel_chain_", chain_name, "_", Sys.Date(), ".rds"))
    } 
  }
  if(verbose == TRUE) close(pb)
  
  #Set colnames for the sample
  if(model$type == "var") {
    if(model$args$sigma == TRUE) {
      colnames(mat) <- c(paste0("A_", 1:(ncol(model$y)*ncol(model$args$xy$xx))),
                         paste0("P_", 1:((ncol(model$y)*(ncol(model$y)+1))/2))
      )
    } else {
      colnames(mat) <- c(paste0("A_", 1:(ncol(model$y)*ncol(model$args$xy$xx))))
    }
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
                 tune = tune,
                 itermax = itermax,
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
  
  model$output <- output
  model
}

