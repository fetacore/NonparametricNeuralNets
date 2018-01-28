#' The Nonparametric polynomial regression with the Least Squares Cross Validation optimal bandwidth
#'
#' @param y   The regressant
#' @param x   The regressors
#' @param h   The bandwidth
#' @param method    The desired kernel function
#' @param p   Order of polynomial local fit
#' @return The vector of fitted values produced by the Leave one out CV
#' @export

nonparpolyLSCV <- function(y, x, h=0.5, method='Gaussian', p=0){

  require(MASS)
  n <- NROW(x)
  g <- matrix(data=NA, nrow=n, ncol=1)
  order = p + 1
  B <- matrix(data=NA, nrow=n, ncol=order)
  e <- matrix(data=NA, nrow=order, ncol=1)
  W <- matrix(data=0, nrow=n, ncol=n)
  k <- kernmatrix(x, h, method)
  LOOy <- matrix(data=0,nrow=n,ncol=1)
  for (i in 1:n){
    for (o in 1:order){
      if (isTRUE(o == 1)){
        e[o] <- 1
      } else {
        e[o] <- 0
      }
      for (j in 1:n){

        if (isTRUE(o-1 == 0)){
          B[j, o] <- 1
        } else {
          B[j, o] <- (x[j]-x[i])^(o - 1)
        }
      }
    }
    LOOy[i] <- y[i]
    diag(W) <- k[i,]
    g[i] <- t(e) %*% ginv(t(B) %*% W %*% B) %*% t(B) %*% W %*% (y-LOOy)
    LOOy[i] <- 0
  }
  LSCV <- 1/n * distance_metric(y,g,'L',2)
  return(list(g,LSCV))
}
