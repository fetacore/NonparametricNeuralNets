#' The general Non parametric regression estimate
#'
#' @param y   The regressant
#' @param x   The observations
#' @param h   The bandwidth
#' @param p     The order of the local polynomial (by default p=0 Nadaraya Watson)
#' @param kernel   The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @return The estimated f for every observation
#' @export



nonparpoly <- function(y, x, h=0.5, p=0, kernel='Gaussian'){
  require(MASS)
  T <- NROW(x)
  f <- matrix(data=NA, nrow=T, ncol=1)
  k <- kernmatrix(x, h, kernel)
  order = p + 1
  B <- matrix(data=NA, nrow=T, ncol=order)
  e <- matrix(data=NA, nrow=order, ncol=1)
  W <- matrix(data=0, nrow=T, ncol=T)
  for (i in 1:T){
    for (o in 1:order){
      if (isTRUE(o == 1)){
        e[o] <- 1
      } else {
        e[o] <- 0
      }
      for (j in 1:T){
        if (isTRUE(o-1 == 0)){
          B[j, o] <- 1
        } else {
          B[j, o] <- (x[j]-x[i])^(o - 1)
        }
      }
    }
    diag(W) <- k[i,]
    f[i] <- t(e) %*% ginv(t(B) %*% W %*% B) %*% t(B) %*% W %*% y
  }
  return(f)
}
