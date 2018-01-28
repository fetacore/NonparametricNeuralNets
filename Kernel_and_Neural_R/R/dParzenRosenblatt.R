#' The first and second derivative of a density based on the Parzen-Rosenblatt density estimator
#'
#' @param x   The observations
#' @param h   The bandwidth
#' @param d   The first or second derivative
#' @note By default we use the Gaussian Kernel
#' @return The estimated f^(r)(x) for every observation
#' @export




dparzen_rosenblatt <- function(x,h,d=1){

  kernel <- 'Gaussian'
  n <- NROW(x)
  k <- diag(x=choose_kernel(0,kernel), nrow=n, ncol=n)
  f <- matrix(data=NA, nrow=n, ncol=1)
  ones <- matrix(data = 1, nrow=n, ncol=1)

  if (isTRUE(d==1)){
    for (i in 1:n){
      for (j in 1:n){
        if (isTRUE(i==j)){
          next
        } else {
          if (isTRUE(k[i,j]==0)){
            k[i,j] <- (-1) * (1/(n*(h)^(2))) * x[i]* choose_kernel(((x[i]-x[j])/h),kernel)
            k[j,i] <- k[i,j]
          } else {
            next
          }
        }
      }
    }
  } else {
    for (i in 1:n){
      for (j in 1:n){
        if (isTRUE(i==j)){
          next
        } else {
          if (isTRUE(k[i,j]==0)){
            k[i,j] <- (1/(n*(h)^(2))) * ((x[i])^2-1)* choose_kernel(((x[i]-x[j])/h),kernel)
            k[j,i] <- k[i,j]
          } else {
            next
          }
        }
      }
    }
  }

  f <- k %*% ones
  return(f)
}
