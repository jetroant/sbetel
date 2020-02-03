
#Generates impulse response functions of a sample from P(A,P,B|Y)
#ADD DOCUMENTATION
irf <- function(model, horizon, narrative = FALSE) {
  
  APB_post <- model$output$APB_sample
  if(narrative == TRUE) {
    APB_post <- model$output$APB_sample_narrative$newSample
  }
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
          irfs[,h,row_index] <- (expm::`%^%`(AA, (h-1)) %*% zero_long)[1:length(zero)]
        }
      }
      
      setTxtProgressBar(pb, row_index)
    }
    close(pb)
    
    ret[[shock_index]] <- irfs
  }
  
  if(narrative == TRUE) {
    model$irfs_narrative <- ret
  } else {
    model$irfs <- ret
  }
  model
}

#Plots impulse response functions
irf_plot <- function(model, narrative = FALSE) {
  
  shocks <- ncol(model$y)
  irfs <- model$irfs
  if(narrative == TRUE) {
    irfs <- model$irfs_narrative
  }
  varnames <- colnames(model$y)
  if(is.null(varnames)) varnames <- paste0("Var. ", 1:shocks)
  
  par(mfrow = c(shocks, shocks))
  indexmat <- matrix(1:shocks^2, ncol = shocks)
  row <- 0
  col <- 0
  for(fig_index in 1:shocks^2) {
    
    if((fig_index-1) %% shocks == 0) {
      row <- row + 1
      col <- 1
    } else {
      col <- col + 1
    }
    
    sub_irfs <- t(irfs[[col]][row,,])
    mean_sub_irfs <- ts(apply(sub_irfs, 2, mean), start = 0)
    
    p <- c(0.049, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.951)
    quant <- function(column) quantile(column, probs = p)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)
    
    color <- "tomato"
    plot(mean_sub_irfs, lwd = 2, lty = 2, col = color, ylab = "", xlab = "", 
         main = paste0("Shock ", col, " on ", varnames[row], " (Nar.)"), ylim = c(min(quantiles_sub_irfs), max(quantiles_sub_irfs)))
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values", probs = p,
        start = 0, fan.col = colorRampPalette(c(color, "white")),
        rlab = NULL, ln = NULL)
    abline(h = 0, lwd = 2, lty = 2)
    if(col == 1 & row == 1) legend("topright", c("90% of post. probability mass"), lwd = 0, bty = "n", col = "tomato")
  }
  par(mfrow = c(1,1))
}