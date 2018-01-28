#' The general Non parametric regression estimate
#'
#' @param y   The regressant
#' @param x   The observations
#' @param h   The bandwidth
#' @param p   The order of the polynomial (either 0 or 1)
#' @param kernel   The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @return The estimated f for every observation
#' @export



multinonparpoly <- function(y, x, h='Silverman', p=1, kernel='Gaussian'){
  require(MASS)
  T <- NROW(x)
  N <- NCOL(x)
  f <- matrix(data=NA, nrow=T, ncol=1)
  h <- the_bandwidth(x, h)
  k <- kernmatrix(x, h, kernel)
  W <- matrix(data=0, nrow=T, ncol=T)
  for (i in 1:T){

    v <- matrix(data=x[i,], nrow=1, ncol=N)
    for (s in 1:N){
      for (j in 1:T){
        if (isTRUE(p==0)){
          B <- matrix(data=1, nrow=T, ncol=1)
        } else {
          B <- matrix(data=1, nrow=T, ncol=N+1)
          B[j,s+1] <- x[j,s] - v[1,s]
        }
      }
    }
    diag(W) <- k[i,]
    f[i] <- sum(ginv(t(B) %*% W %*% B) %*% t(B) %*% W %*% y)
  }
  return(f)
}
