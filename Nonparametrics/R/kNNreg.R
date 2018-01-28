#' The k-nearest neighborhood regression
#'
#' @param y The regressant
#' @param x The regressor
#' @param k The number of neighbors (1 neighbor by default)
#' @param weight Assigns different weights to neighbors according to distance
#' @param method If weight=TRUE then this is the kernel function to evaluate the weights
#' @export

kNNreg <- function(y,x,k=1,weight=FALSE, method='Gaussian'){

  T <- NROW(x)
  a <- matrix(data=NA, nrow=T, ncol=1)
  b <- matrix(data=NA, nrow=T, ncol=T)
  f <- matrix(data=NA, nrow=T, ncol=1)
  s <- matrix(data=NA, nrow=k, ncol=1)
  kern <- matrix(data=NA, nrow=k, ncol=1)
  for (i in 1:T){
    a[i] <- i
    for (j in 1:T){
      b[j,i] <- abs(x[i]-x[j])
    }
  }
  d <- matrix(c(a, b), nrow=T, ncol=T+1)
  if (isTRUE(weight)){
    for (i in 1:T){
      o <- matrix(c(d[,1],d[,i+1]), nrow=T, ncol=2)
      o <- o[order(o[,2]),]
      for (j in 1:k){
        if (isTRUE(j==1)){
          kern[j] <- choose_kernel(0, method)
          s[j] <- choose_kernel(0, method) * y[o[1,1]]
        } else {
          kern[j] <- choose_kernel((x[i]-x[o[j,1]])/which.min(o[2:NROW(o),2]),method)
          s[j] <- kern[j] * y[o[j,1]]
        }
      }
      f[i] <- sum(s)/sum(kern)
    }
  } else {
    for (i in 1:T){
      if (!isTRUE(k==T)){
        o <- matrix(c(d[,1],d[,i+1]), nrow=T, ncol=2)
        o <- o[order(o[,2]),]
        f[i] <- 1/k * sum(y[o[1,1]:o[k+1,1]])
      } else {
        o <- matrix(c(d[,1],d[,i+1]), nrow=T, ncol=2)
        o <- o[order(o[,2]),]
        f[i] <- 1/T * sum(y[o[1,1]:o[k,1]])
      }
    }
  }
  return(f)
}
