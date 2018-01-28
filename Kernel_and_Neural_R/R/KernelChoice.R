#'	Defines the desired kernel for our problem.
#'
#'	@param u	A	number/vector.
#'	@param arg The optional kernel (Alphabetical)
#'	@note If argument is not alphabetical, it returns the Gaussian Kernel
#'	@importFrom stats cov
#'	@return The	kernelized values.
#'	@export


choose_kernel <- function(u, arg=NULL){

  k <- rep(0, NROW(u))

  if (isTRUE(arg == 'Epanechnikov')){

    for (i in 1:NROW(u)){
      if (abs(u[i]) <= 1){
        I <- 1
      }else {
        I <- 0
      }
      k[i] <- (3/4)*(1 - (u[i])^2)*I
    }
    kernel <- k
  } else if (isTRUE(arg == 'Silverman')){

    for (i in 1:NROW(u)){
      k[i] <- (1/2)*exp(-(abs(u[i])/sqrt(2)))*sin((abs(u[i])/sqrt(2))+(pi/4))
    }
    kernel <- k
  } else if (isTRUE(arg == 'Biweight')){

    for (i in 1:NROW(u)){
      if (abs(u[i]) <= 1){
        I <- 1
      }else {
        I <- 0
      }
      k[i] <- (15/16)* (1 - (u[i])^2)^2 * I
    }
    kernel <- k
  } else {

    for (i in 1:NROW(u)){
      k[i] <- (1/sqrt(2 * pi))*exp(-(u[i]^2) / 2)
    }
    kernel <- k
  }
  return(kernel)
}
