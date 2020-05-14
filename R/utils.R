
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



