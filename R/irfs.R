
#Generates impulse response functions of a sample from P(A,P,B|Y)
#ADD DOCUMENTATION
irf <- function(APB_post, horizon) {
    
  A_post <- APB_post[,grep("A", colnames(APB_post))]
  B_post <- APB_post[,grep("B", colnames(APB_post))]
  if(sqrt(ncol(B_post)) != floor(sqrt(ncol(B_post)))) stop("Something is wrong with 'APB_post'")
  shocks <- sqrt(ncol(B_post))
  
  ret <- list()
  
  for(shock_index in 1:shocks) {
    
    e <- rep(0, shocks)
    e[shock_index] <- 1
    irfs <- array(NA, dim = c(shocks, horizon+1, nrow(A_post)))
    
    print(paste0("Computing impulse responses... (", shock_index, "/", shocks,")"))
    pb <- txtProgressBar(min = 0, max = nrow(A_post), style = 3)
    for(row_index in 1:nrow(A_post)) {
      
      B <- matrix(B_post[row_index,], ncol = shocks)
      A <- matrix(A_post[row_index,], ncol = shocks)
      AA <- stackA(A)
      
      for(h in 1:(horizon+1)) {
        if(h == 1) {
          zero <- B %*% e
          zero_long <- c(zero, rep(0, (ncol(AA)-length(zero))))
          irfs[,h,row_index] <- zero
        } else {
          irfs[,h,row_index] <- (AA %^% (h-1) %*% zero_long)[1:length(zero)]
        }
      }
      
      setTxtProgressBar(pb, row_index)
    }
    close(pb)
    
    ret[[shock_index]] <- irfs
  }
  ret
}
