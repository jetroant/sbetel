
#Builds xx and yy from y and p (VAR)
build_xy_var <- function(y, p) {
  n <- ncol(y)
  for(i in 1:p) {
    if(i == 1) {
      temp <- rbind(rep(NA, n), y[-(nrow(y)),])
      xx <- temp
    } else {
      temp <- rbind(rep(NA, n), temp[-(nrow(y)),])
      xx <- cbind(xx, temp)
    }
  }
  xx <- cbind(rep(1, nrow(xx)), xx)
  xx <- xx[-c(1:p),]
  yy <- y[-c(1:p),]
  ret <- list("xx" = xx, "yy" = yy)
  ret
}

#Returns all permutations of 'n' integers
permutations <- function(n){
  if(n == 1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow = n*p, ncol = n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i, sp + (sp >= i))
    }
    return(A)
  }
}

#Load backup folder
load_backup <- function(dir) {
  
  #Load the backup directory
  wd <- getwd()
  setwd(dir)
  subchains <- list()
  filenames <- list.files()
  for(i in 1:length(filenames)) {
    subchains[[i]] <- readRDS(filenames[i])
  }
  setwd(wd)
  
  #Set parameters
  chain_length <- nrow(subchains[[1]]$chain)
  chain_number <- length(subchains)
  
  #Collect the chains and likelihoods
  allchains <- matrix(NA, 
                      ncol = ncol(subchains[[1]]$chain),
                      nrow = chain_length*chain_number
  )
  likelihoods <- c()
  for(i in 1:length(subchains)) {
    rows <- (i*chain_length+1):((i+1)*chain_length)-chain_length
    allchains[rows,] <- subchains[[i]]$chain
    burn <- length(subchains[[i]]$likelihoods) - nrow(subchains[[i]]$chain)
    if(burn > 0) {
      likelihoods <- c(likelihoods, subchains[[i]]$likelihoods[-c(1:burn)])
    } else {
      likelihoods <- c(likelihoods, subchains[[i]]$likelihoods)
    }
  }
  
  #Average acceptance rate
  accrates <- rep(NA, chain_number)
  for(i in 1:chain_number) accrates[i] <- subchains[[i]]$accrate
  
  #Collect the output
  ret <- list(sample = allchains,
              subchains = subchains,
              accrates = accrates,
              likelihoods = likelihoods,
              totaltime = NA,
              tune = NA,
              itermax = NA,
              type = NA)
  ret
}


