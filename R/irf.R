
irf <- function(model, output, 
                horizon = 40, N = 10000,
                cumulate = c()) {
  
  m <- ncol(model$y)
  p <- model$args$p
  
  A_indices <- 1:(m*(1+m*p))
  B_indices <- (m*(1+m*p)+1):(m*(1+m*p)+m^2) 
  
  ret <- list()
  rows <- sample.int(nrow(output$sample), N, replace = TRUE)
  A_post <- output$sample[rows, A_indices]
  B_post <- output$sample[rows, B_indices]
  
  for(shock_index in 1:m) {
    
    e <- rep(0, m)
    e[shock_index] <- 1
    irfs <- array(NA, dim = c(m, horizon+1, nrow(A_post)))
    
    print(paste0("Computing impulse responses... (", shock_index, "/", m,")"))
    pb <- txtProgressBar(min = 0, max = nrow(A_post), style = 3)
    for(row_index in 1:nrow(A_post)) {
      
      B <- matrix(B_post[row_index,], ncol = m)
      A <- matrix(A_post[row_index,], ncol = m)
      AA <- sbetel:::stackA(A)
      
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
    
    if(length(cumulate) > 0) {
      for(i in 1:length(cumulate)) {
        for(j in 1:N) {
          irfs[cumulate[i],,j] <- cumsum(irfs[cumulate[i],,j])
        }
      }
    }
    
    ret[[shock_index]] <- irfs
  }
  ret
}

#Stacks VAR(p) coefficient matrix to VAR(1) coefficient matrix
stackA <- function(A) {
  A <- t(A)[,-1]
  m <- nrow(A)
  lags <- ncol(A)/m
  eye <- diag(m*lags-m)
  A <- rbind(A, cbind(eye, matrix(0, ncol = m, nrow= nrow(eye))))
  A
}

#Plots impulse response functions
irf_plot <- function(irf_obj, varnames, probs = NULL) {
  
  m <- nrow(irf_obj[[1]]) 
  par(mar = c(2,4,2,1))
  par(mfrow = c(m, m))
  indexmat <- matrix(1:m^2, ncol = m)
  row <- 0
  col <- 0
  for(fig_index in 1:m^2) {
    
    if((fig_index-1) %% m == 0) {
      row <- row + 1
      col <- 1
    } else {
      col <- col + 1
    }
    
    sub_irfs <- t(irf_obj[[row]][col,,])
    mean_sub_irfs <- ts(apply(sub_irfs, 2, mean), start = 0)
    
    if(is.null(probs)) {
      p <- c(0.0249, 0.025, seq(0.1, 0.9, 0.1), 0.975, 0.9751)
    } else {
      p <- probs
    }
    quant <- function(column) quantile(column, probs = p)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)
    
    color <- "tomato"
    plot(mean_sub_irfs, lwd = 2, lty = 2, col = color, ylab = "", xlab = "", 
         main = paste0("Shock ", row, " on ", varnames[col]), 
         ylim = c(min(quantiles_sub_irfs), max(quantiles_sub_irfs)))
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values", probs = p,
                 start = 0, fan.col = colorRampPalette(c(color, "white")),
                 rlab = NULL, ln = NULL)
    abline(h = 0, lwd = 2, lty = 2)
    
    post_mass <- (max(p[-length(p)]) - min(p[-1]))*100
    if(col == 1 & row == 1) legend("topright", c(paste0(post_mass,"% of post. prob. mass")), lwd = 0, bty = "n", col = "tomato")
  }
  par(mfrow = c(1,1))
}


