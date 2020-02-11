
#Draws 'PQ'/'B' matrix that satisfies the given sign restrictions in 'res'
PQ_try <- function(P, res, perms = permutations(ncol(P)), additional = FALSE) {
  
  k <- ncol(P)
  W <- matrix(rnorm(k^2, 0, 1), ncol = k)
  qr_obj <- qr(W)
  Q <- qr.Q(qr_obj) %*% diag(ifelse(diag(qr.R(qr_obj)) < 0, -1, 1))
  
  accept <- FALSE
  
  PQ <- P %*% Q
  PQ_signs <- PQ / abs(PQ)
  
  #Signs match
  if(sum(PQ_signs != res) == 0) {
    accept <- TRUE
    #print("first")
    
    #Opposite signs match  
  } else if(sum(PQ_signs == res) == 0) {
    PQ <- P %*% (Q * (-1))
    accept <- TRUE
    #print("opposite")
    
    #Permutations of Q
  } else {
    
    for(i in 2:nrow(perms)) {
      Q_perm <- Q[,perms[i,]]
      PQ <- P %*% Q_perm
      PQ_signs <- PQ / abs(PQ)
      
      #Signs match
      if(sum(PQ_signs != res) == 0) {
        accept <- TRUE
        #print("perm first")
        
        #Opposite signs match  
      } else if(sum(PQ_signs == res) == 0) {
        PQ <- P %*% (Q_perm * (-1))
        accept <- TRUE
        #print("perm opposite")
      }
      
      #If accepted, stop the permutation search
      if(accept == TRUE) break
    }
    
  }
  
  #Additional sign restrictions
  if(additional == TRUE) {
    
    #Price elasticity of oil supply
    if(PQ[1,2]/PQ[3,2] > 0.0258 | PQ[1,3]/PQ[3,3] > 0.0258) {
      accept <- FALSE
    }
  }
  
  if(accept == TRUE) {
    return(PQ)
  } else {
    return(NULL)
  }
}

#Generates a sample (of size 'N') from the posterior P(A,P,B|Y) given 
#a sample from P(A,P|Y) and the sign restrictions in 'res'
#(ADD DOCUMENTATION)
B_signres_sample <- function(model, 
                             res, 
                             N, 
                             additional = FALSE, 
                             gaussian = FALSE,
                             verbose = TRUE,
                             backup = NULL,
                             core_number = 0) {
  
  shocks <- ncol(model$y)
  AP_post <- model$output$sample
  if(gaussian == TRUE) AP_post <- model$output$gaussian_sample
  A_post <- model$output$sample[,grep("A", colnames(AP_post))]
  P_post <- model$output$sample[,grep("P", colnames(AP_post))]
  
  APB_post <- matrix(NA, ncol = ncol(A_post) + ncol(P_post) + shocks^2, nrow = N)
  colnames(APB_post) <- c(paste0("A_", 1:ncol(A_post)), 
                          paste0("P_", 1:ncol(P_post)), 
                          paste0("B_", 1:(shocks^2)))
  
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    
    PQ <- NULL
    while(is.null(PQ)) {
      P <- diag(shocks)
      post_index <- sample.int(nrow(P_post),1)
      P[!upper.tri(P)] <- P_post[post_index,]
      PQ <- PQ_try(P, res, additional = additional)
    }
    
    APB_post[i,1:ncol(A_post)] <- A_post[post_index,]
    APB_post[i,(ncol(A_post)+1):(ncol(A_post)+ncol(P_post))] <- P_post[post_index,]
    APB_post[i,(ncol(A_post)+ncol(P_post)+1):ncol(APB_post)] <- c(PQ)
    if(verbose == TRUE) setTxtProgressBar(pb, i)
    
    if(!is.null(backup)) {
      if(i %% backup == 0) saveRDS(list(mat, likelihoods), paste0("sbetel_signres_sample_", core_number, "_", Sys.Date(), ".rds"))
    } 
  }
  if(verbose == TRUE) close(pb)
  
  if(gaussian == TRUE) model$output$APB_sample_gaussian <- APB_post
  if(gaussian == FALSE) model$output$APB_sample <- APB_post
  model
}


