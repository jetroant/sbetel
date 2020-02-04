
#ADD DOCUMENTATION
est_sbetel_parallel <- function(model,
                                chains = parallel::detectCores() - 1,
                                N,
                                tune = "auto",
                                itermax = 20, 
                                burn = NULL,
                                n = 100,
                                verbose = TRUE,
                                backup = 100) {
  
  #Choosing the tuning parameter for suitable step size
  if(tune == "auto") {
    autotune_output <- auto_tune(model, n = n, verbose = verbose)
    tune <- autotune_output$tune
  }
  
  #Initiate cluster
  cl <- parallel::makeCluster(chains)
  
  #Export the needed objects/arguments
  parallel::clusterExport(cl, list("model", "N", "tune"))
  
  #Run parallel chains
  model_copies <- parallel::parLapply(cl,
                                      1:chains,
                                      function(chain_number) {
                                        sbetel::est_sbetel(model = model,
                                                           N = N,
                                                           tune = tune,
                                                           backup = 1000,
                                                           chain_name = chain_number,
                                                           verbose = FALSE)
                                      })
  
  parallel::stopCluster(cl)
  
  #Collect the chains
  model$output_parallel <- model_copies
  
  #Re-save and overwrite old sample (if any)
  if(!is.null(model$output$sample)) {
    model$output$sample_old <- model$output$sample
  }
  for(i in 1:chains) {
    if(i == 1) {
      sample_all <- model_copies[[i]]$output$sample[-1,]
    } else {
      sample_all <- rbind(sample_all, model_copies[[i]]$output$sample[-1,])
    }
  }
  model$output$sample <- sample_all
  
  model
}



