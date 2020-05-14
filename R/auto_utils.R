
#Approximates efficient step size for RWMH algortihm,
#mimicking the manual tuning process.
#Gets you in the right ballpark.
auto_tune <- function(model,
                      type,
                      target = c(0.22, 0.28), 
                      limit = 20, 
                      chain_length = 100, 
                      verbose = TRUE) {
  
  if(verbose == TRUE) cat(paste0("Auto tuning step size... \n"))
  found <- FALSE
  tune <- 2.38^2/length(model$initial(model$args$xy)$th)
  step <- 0.1
  direction <- c()
  while(found == FALSE) {
    
    #Autotune fails if limit is reached
    if(length(direction) > limit) stop("Autotune failed. Limit reached.")
    
    s <- est_sbetel(model,
                    chain_length = chain_length,
                    chain_number = 1,
                    tune = tune,
                    burn = 0,
                    parallel = 1,
                    verbose = FALSE)
    if(verbose == TRUE) cat(paste0("tune = ", round(tune, 4), " / accrate = ", round(s$accrates, 2), "\n"))
    
    if(s$accrates > target[2]) {
      
      #Step size halves if the direction of the search changes
      direction <- c(1, direction)
      if(length(direction) != 1 & direction[1] == -direction[2]) step <- step/2
      tune <- tune + step
    } else if(s$accrates < target[1]) {
      
      #Step size halves if the direction of the search changes
      direction <- c(-1, direction)
      if(length(direction) != 1 & direction[1] == -direction[2]) step <- step/2
      if((tune - step) < 0) {
        tune <- tune/2
      } else {
        tune <- tune - step
      }
      
    } else {
      
      #Parameter is accepted if 'accrate' is still within the target
      s <- est_sbetel(model,
                      chain_length = chain_length,
                      chain_number = 1,
                      tune = tune,
                      burn = 0,
                      parallel = 1,
                      verbose = FALSE)
      if(s$accrates > target[1] & s$accrates < target[2]) {
        found <- TRUE
      }
    }
  }
  
  list(accrate = s$accrates, tune = tune, iters = length(direction))
}

