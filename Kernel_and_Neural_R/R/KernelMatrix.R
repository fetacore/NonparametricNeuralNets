#' The Kernel Matrix
#'
#' @param x   The observations
#' @param h  The bandwidth
#' @param kernel The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @note By default we use the Gaussian Kernel
#' @return The k((x[i]-x[j])/h) matrix with i rows and j columns
#' @export




kernel_matrix <- function(x, h, kernel='Gaussian'){

  T <- NROW(x)
  k <- matrix(data=NA, nrow=T, ncol=T)
  for (i in 1:T){
    for (j in 1:T){
      k[i,j] <- choose_kernel(((x[i]-x[j])/h),kernel)
    }
  }
  return(k)
}
