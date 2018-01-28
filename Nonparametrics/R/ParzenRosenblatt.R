#' The Parzen-Rosenblatt density estimator
#'
#' @param x   The observations
#' @param h   The bandwidth
#' @param kernel The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @note By default we use the Gaussian Kernel
#' @return The estimated f for every observation
#' @export




parzen_rosenblatt <- function(x,h,kernel='Gaussian'){

  T <- NROW(x)
  k <- diag(x=choose_kernel(0,kernel), nrow=T, ncol=T)
  f <- matrix(data=NA, nrow=T, ncol=1)
  ones <- matrix(data = 1, nrow=T, ncol=1)

  for (i in 1:T){
    for (j in 1:T){
      if (isTRUE(i==j)){
        next
      } else {
        if (isTRUE(k[i,j]==0)){
          k[i,j] <- (1/(T*h)) * choose_kernel(((x[i]-x[j])/h),kernel)
          k[j,i] <- k[i,j]
        } else {
          next
        }
      }
    }
  }
  f <- k %*% ones
  return(f)
}
