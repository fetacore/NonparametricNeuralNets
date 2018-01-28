#' The general Kernel Matrix
#'
#' @param x   The observations
#' @param h  The bandwidth
#' @param distance The distance metric
#' @param kernel The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @note By default we use the Gaussian Kernel
#' @return The k((x[i]-x[j])/h) matrix with i rows and j columns
#' @export




gkernel_matrix <- function(x, h, distance='Eucledean', kernel='Gaussian'){

  T <- NROW(x)
  N <- NCOL(x)
  k <- matrix(data=NA, nrow=T, ncol=T)
  for (i in 1:T){
    if (isTRUE(N==1)){
      v <- x[i]
      for (j in 1:T){
        k[i,j] <- choose_kernel((distance_metric(v,x[j],distance)/h),kernel)
      }
    } else {
      v <- matrix(data=x[i,], nrow=1, ncol=N)
      for (j in 1:T){
        k[i,j] <- choose_kernel((distance_metric(v,x[j,],distance)/h),kernel)
      }
    }
  }
  return(k)
}
