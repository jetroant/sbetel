
narratives <- function(model, 
                       output,
                       res_fun, 
                       N = 1000,
                       M = 10000,
                       start_date = NULL,
                       freq = NULL,
                       parallel = 1,
                       verbose = TRUE) {
  
  y <- model$y
  yy <- model$args$xy$yy
  xx <- model$args$xy$xx
  if(!("ts" %in% class(y))) { 
    if(is.null(start_date) | is.null(freq)) {
      stop("'start_date' and 'freq' need to provided if 'model$y' is not a time series object.")
    } 
    y <- ts(y, start = start_date, frequency = freq)
  } else {
    freq <- attributes(y)$tsp[3]
  }
  start_date_e <- attributes(y)$tsp[1] + model$args$p/freq 
  
  A_cols <- 1:(ncol(y)*(1+ncol(y)*model$args$p))
  B_cols <- (ncol(y)*(1+ncol(y)*model$args$p)+1):(ncol(y)*(1+ncol(y)*model$args$p)+ncol(y)^2)
  
  s <- output$sample
  accepted <- rep(FALSE, nrow(s))
  if(verbose == TRUE) cat("Checking narrative restrictions... (1/2) \n")
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = nrow(s), style = 3) 
  for(i in 1:nrow(s)) {
    
    A <- matrix(s[i,A_cols], ncol = ncol(y))
    B <- matrix(s[i,B_cols], ncol = ncol(y))
    B_inv <- solve(B)
    U <- yy - xx %*% A
    E <- ts(U %*% t(B_inv), start = start_date_e, frequency = freq)
    
    accepted[i] <- res_fun(E = E, A = A, B = B)
    if(verbose == TRUE) setTxtProgressBar(pb, i)
  }
  if(verbose == TRUE) close(pb)
  acc_rate <- (sum(accepted)/length(accepted))*100
  if(verbose == TRUE) cat(paste0(round(acc_rate, 3), "% of the posterior sample accepted. \n", 
                                 "(", sum(accepted), " / ", length(accepted), " draws) \n"))
  
  if(sum(accepted) == 0) {
    if(verbose == TRUE) cat("Narrative restrictions too tight, NULL returned.")
    return(NULL)
  }
  
  s <- s[which(accepted == TRUE),]
  importance <- rep(NA, nrow(s))
  if(verbose == TRUE) cat("Calculating importance weights... (2/2) \n")
  if(parallel == 1) {
    if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = nrow(s), style = 3) 
    for(i in 1:nrow(s)) {
      importance[i] <- sbetel:::calculate_importance(s_row = s[i,], res_fun, yy, xx, 
                                                     A_cols, B_cols, start_date_e, freq, M)
      if(verbose == TRUE) setTxtProgressBar(pb, i)
    }
    if(verbose == TRUE) close(pb)
  }
  if(parallel > 1) {
    cat(paste0("Calculating in parallel, using ", parallel, " cores...\n"))
    cl <- parallel::makeCluster(parallel)
    parallel::clusterExport(cl,
                            list("s", "yy", "xx", "A_cols", "B_cols",
                                 "start_date_e", "freq", "M", "res_fun"),
                            envir = environment())
    importance <- parallel::parSapply(cl,
                                      1:nrow(s),
                                      function(i) {
                                        sbetel:::calculate_importance(s_row = s[i,],
                                                                      res_fun,
                                                                      yy, 
                                                                      xx, 
                                                                      A_cols,
                                                                      B_cols,
                                                                      start_date_e,
                                                                      freq,
                                                                      M)
                                      })
    parallel::stopCluster(cl)
  }
  importance[which(importance == Inf)] <- M
  
  #Resampling according to importance weights
  s_new <- s[sample.int(nrow(s), N, replace = TRUE, prob = importance),]
  
  ret <- list("sample" = s_new,
              "raw_sample" = s,
              "importance_weights" = importance)
  ret
}

calculate_importance <- function(s_row, res_fun, yy, xx, 
                                 A_cols, B_cols, start_date_e, freq, M) {
  
  A <- matrix(s_row[A_cols], ncol = ncol(yy))
  B <- matrix(s_row[B_cols], ncol = ncol(yy))
  B_inv <- solve(B)
  U <- yy - xx %*% A
  E <- ts(U %*% t(B_inv), start = start_date_e, frequency = freq)
  
  sub_accepted <- rep(NA, M)
  for(j in 1:M) { 
    E_rndm <- ts(E[sample.int(nrow(E), nrow(E), replace = T),],
                 start = start_date_e, frequency = freq)
    sub_accepted[j] <- res_fun(E = E_rndm, A = A, B = B)
  }
  
  1/mean(sub_accepted)
}



