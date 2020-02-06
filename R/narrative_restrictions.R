
#util for 'narrative_sample()' 1/5
calculateH <- function(Ls, E, periods, variable, yy) {
  H <- matrix(NA, ncol = ncol(yy), nrow = length(periods))
  for(period in periods) {
    h <- period - min(periods)
    
    #Diaz & Rubio-Ramirez, 2018, p.2805
    for(shock in 1:ncol(yy)) {
      H[(period-min(periods))+1, shock] <- t(diag(ncol(yy))[,variable]) %*% Ls[[h+1]] %*% diag(ncol(yy))[,shock] %*% t(diag(ncol(yy))[,shock]) %*% E[max(periods)-h,]
      
    }
  }
  apply(H, 2, sum)
}

#util for 'narrative_sample()' 2/5
histDecomp <- function(variable, periods, param, xx, yy) {
  
  horizon <- max(periods) - min(periods)
  
  A <- matrix(param[grep("A", names(param))], ncol = ncol(yy))
  B <- matrix(param[grep("B", names(param))], ncol = ncol(yy))
  U <- yy - xx %*% A
  E <- U %*% t(solve(B))
  AA <- stackA(A)
  
  if(variable != 0) {
    
    #Initialize Ls
    Ls <- list()
    for(i in 1:(horizon+1)) {
      Ls[[i]] <- matrix(NA, ncol = ncol(yy), nrow = ncol(yy))
    }
    
    #Calculate Ls
    for(j in 1:ncol(yy)) {
      e <- rep(0, ncol(yy))
      e[j] <- 1
      for(h in 1:(horizon+1)) {
        if(h == 1) {
          zero <- B %*% e
          zero_long <- c(zero, rep(0, (ncol(AA)-length(zero))))
          Ls[[h]][,j] <- zero
        } else {
          Ls[[h]][,j] <- (expm::`%^%`(AA, (h-1)) %*% zero_long)[1:length(zero)]
        }
      }
    }
    
    #Calculate H
    H <- calculateH(Ls, E, periods, variable, yy)
    ret <- list(H = H, Ls = Ls, signs = NULL, periods = periods, variable = variable)
  } else {
    
    #Check signs
    signs <- E[periods,]
    signs <- matrix(signs/abs(signs), nrow = length(periods))
    ret <- list(H = NULL, Ls = NULL, signs = signs, periods = periods, variable = variable)
  }
  
  ret
}

#util for 'narrative_sample()' 3/5
eFun <- function(param, xx, yy) {
  A <- matrix(param[grep("A", names(param))], ncol = ncol(yy))
  B <- matrix(param[grep("B", names(param))], ncol = ncol(yy))
  U <- yy - xx %*% A
  E <- U %*% t(solve(B))
  E
}

#util for 'narrative_sample()' 4/5
H_accept <- function(H_obj, narrative) {
  
  if(sum(narrative != 0) != 1) stop("Narrative not properly defined.")
  H <- H_obj$H
  signs <- H_obj$signs[,which(narrative != 0)]
  accept <- FALSE
  
  #"... shock most important contributor"
  if(sum(narrative) == 1) {
    if(abs(H[which(narrative != 0)]) > max(abs(H[which(narrative == 0)]))) accept <- TRUE
    
    #"... shock overwhelming contributor
  } else if(sum(narrative) == 2) {
    if(abs(H[which(narrative != 0)]) > sum(abs(H[which(narrative == 0)]))) accept <- TRUE
    
    #"... shock least important contributor
  } else if(sum(narrative) == -1) {
    if(abs(H[which(narrative != 0)]) < min(abs(H[which(narrative == 0)]))) accept <- TRUE
    
    #... shock NOT overwhelming contributor
  } else if(sum(narrative) == -2) {
    if(abs(H[which(narrative != 0)]) < sum(abs(H[which(narrative == 0)]))) accept <- TRUE
    
    #Sign of shocks positive
  } else if(sum(narrative) == 3) {
    if(sum(signs) == length(signs)) accept <- TRUE
    
    #Sign of shocks negative
  } else if(sum(narrative) == -3) {
    if(sum(signs) == -length(signs)) accept <- TRUE
    
  }
  
  accept
}

#util for 'narrative_sample()' 5/5
importanceWeight <- function(H_objects, narrative_obj, E, xy, n = 1000, p = p, gaussian = FALSE) {
  
  yy <- xy$yy
  td <- ifelse(is.null(xy$td), 0, xy$td)
  
  #Periods for which shocks are resampled
  periods <- c()
  for(i in 1:length(narrative_obj$periods)) periods <- c(periods, narrative_obj$periods[[i]] + td - p)
  periods <- which(c(1:nrow(E)) %in% periods)
  
  E_new <- E
  avec <- rep(TRUE, n)
  for(i in 1:n) {
    
    #Resample the structural shocks (in RR 2018 iid standard normal assumed, i.e. gaussian == TRUE)
    if(gaussian == FALSE) {
      E_new[periods,] <- E[sample((td+1):nrow(E), length(periods), replace = TRUE),]
    } else {
      E_new[periods,] <- matrix(rnorm(length(periods)*ncol(E)), ncol = ncol(E))
    }

    accept <- TRUE
    for(j in 1:length(narrative_obj$restrictions)) {
      
      if(narrative_obj$variables[[j]] != 0) {
        
        #Calculate new H
        H <- calculateH(Ls = H_objects[[j]]$Ls, E = E_new, periods = H_objects[[j]]$periods + td, variable = H_objects[[j]]$variable, yy = yy)
        H_obj_new <- list(H = H, Ls = H_objects[[j]]$Ls, signs = NULL, periods = H_objects[[j]]$periods + td, variable = H_objects[[j]]$variable)
        
        #Check if accepted
        accept <- H_accept(H_obj_new, narrative_obj$restrictions[[j]])
        if(accept == FALSE) break
        
      } else {
        
        #Check signs
        signs <- E_new[H_objects[[j]]$periods + td,]
        signs <- matrix(signs/abs(signs), nrow = length(H_objects[[j]]$periods + td))
        H_obj_new <- list(H = NULL, Ls = NULL, signs = signs, periods = H_objects[[j]]$periods + td, variable = H_objects[[j]]$variable)
        
        #Check if accepted
        accept <- H_accept(H_obj_new, narrative_obj$restrictions[[j]])
        if(accept == FALSE) break
      }
      
    }
    
    avec[i] <- accept
  }
  
  #Importance weight
  1 -mean(avec)
}

#Generates a sample from the posterior P(A,P,B|Y,Nar)
#conditional on narrative restrictions given a sample 
#from P(A,P,B|Y) and the narrative restrictions in 'narrative_obj'
narrative_sample <- function(model, narrative_obj, 
                             N = nrow(model$output$APB_sample),
                             gaussian = FALSE) {
  
  if(gaussian == FALSE) APB_post <- model$output$APB_sample
  if(gaussian == TRUE) APB_post <- model$output$APB_sample_gaussian
  xy <- model$args$xy
  xx <- xy$xx
  yy <- xy$yy
  td <- ifelse(is.null(xy$td), 0, xy$td)
  p <- model$args$p
  
  acceptedSample <- matrix(c(APB_post[1,], NA), nrow = 1)
  acceptedSample[1,] <- NA
  colnames(acceptedSample) <- c(colnames(APB_post), "ImportanceWeight")
  
  #Go through the whole sample
  print("Checking narrative restrictions... (1/2)")
  pb <- txtProgressBar(min = 0, max = nrow(APB_post), style = 3)
  for(i in 1:nrow(APB_post)) {
    
    #Accept or discard a draw 
    accept <- TRUE
    H_objects <- list()
    for(j in 1:length(narrative_obj$restrictions)) {
      
      H_objects[[j]] <- histDecomp(variable = narrative_obj$variables[[j]], 
                                   periods = narrative_obj$periods[[j]] + td - p, 
                                   param = APB_post[i,], 
                                   xx = xx, 
                                   yy = yy)
      accept <- H_accept(H_obj = H_objects[[j]], narrative = narrative_obj$restrictions[[j]])
      if(accept == FALSE) break
    }
    if(accept == TRUE) {
      
      #Calculate the importance weight
      E <- eFun(param = APB_post[i,], xx = xx, yy = yy)
      iw <- importanceWeight(H_objects = H_objects, narrative_obj = narrative_obj, 
                             E = E, xy = xy, n = 1000, p = p, gaussian = gaussian)
      acceptedSample <- rbind(acceptedSample, 
                              c(APB_post[i,], iw))
    }
    setTxtProgressBar(pb, i)
  } 
  close(pb)
  acceptedSample <- acceptedSample[-1,]
  #If only one row this needs to be done in order for the sample not to become a vector
  if(is.null(nrow(acceptedSample))) {
    cnames <- names(acceptedSample)
    acceptedSample <- data.frame(matrix(acceptedSample, nrow = 1))
    colnames(acceptedSample) <- cnames
  }
  
  print(paste0(100*nrow(acceptedSample)/nrow(APB_post), "% of the sample accepted (Narrative restrictions)..."))
  if(nrow(acceptedSample) == 0) {
    print("Narrative restrictions too tight. No new sample returned.")
    return(model)
  }
  
  #Resample using importance weights
  newSample <- matrix(NA, ncol = ncol(APB_post), nrow = N)
  iws <- acceptedSample[,"ImportanceWeight"]
  iws_cumsum <- cumsum(iws)
  unidraw <- runif(N, 0, max(iws_cumsum))
  print("Resampling using importance weights... (2/2)")
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    row_index <- which(iws_cumsum == min(iws_cumsum[which(iws_cumsum >= unidraw[i])]))[1]
    if(nrow(acceptedSample) != 1) the_row <- acceptedSample[row_index, -ncol(acceptedSample)]
    if(nrow(acceptedSample) == 1) the_row <- as.numeric(acceptedSample[row_index, -ncol(acceptedSample)])
    newSample[i,] <- the_row
    setTxtProgressBar(pb, i)
  }
  close(pb)
  print("Done!")
  
  colnames(newSample) <- colnames(APB_post)
  ret <- list(newSample = newSample, acceptedSample = acceptedSample)
  if(gaussian == FALSE) model$output$APB_sample_narrative <- ret
  if(gaussian == TRUE) model$output$APB_sample_gaussian_narrative <- ret
  
  model
}
