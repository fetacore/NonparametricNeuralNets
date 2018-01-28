#' The Kernel Matrix
#'
#' @param x   The observations (vector or matrix)
#' @param h  The bandwidth
#' @param kernel The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @note By default we use the Gaussian Kernel
#' @return The k((x[i]-x[j])/h) matrix with i rows and j columns
#' @export




kernmatrix <- function(x, h, kernel='Gaussian'){

  T <- NROW(x)
  N <- NCOL(x)
  kernmatrix <- matrix(data=NA, nrow=T, ncol=T)
  for (i in 1:T){
    if (isTRUE(N==1)){
      v <- x[i]
      for (j in 1:T){
        kernmatrix[i,j] <- choice_of_kernel(((v-x[j])/h),kernel)
      }
    } else {
      v <- matrix(data=x[i,], nrow=1, ncol=N)
      for (j in 1:T){
        kernmatrix[i,j] <- choice_of_kernel(((v-x[j,])/h),kernel)
      }
    }
  }
  return(kernmatrix)
}
