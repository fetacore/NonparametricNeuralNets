#'	Defines the desired kernel for our problem.
#'
#'	@param x	A	number/row vector.
#'	@param arg The optional kernel (Alphabetical)
#'	@note If argument is not alphabetical, it returns the Gaussian Kernel
#'	@importFrom stats cov
#'	@return The	kernelized values.
#'	@export


choice_of_kernel <- function(x, arg=NULL){

#Row vector argument

  if (isTRUE(x-x == matrix(data=0,nrow=NROW(x),ncol=NCOL(x)))){

    N <- NCOL(x)
    k <- matrix(data=NA, nrow=T,ncol=N)

    if (isTRUE(arg == 'Epanechnikov')){

      for (i in 1:N){
        if (isTRUE(abs(x[i]) <= 1)){
          I <- 1
        }else {
          I <- 0
        }
        k[i] <- (3/4)*(1 - (x[i])^2)*I
      }

    } else if (isTRUE(arg == 'Silverman')){

      for (i in 1:N){
        k[i] <- (1/2)*exp(-(abs(x[i])/sqrt(2)))*sin((abs(x[i])/sqrt(2))+(pi/4))
      }

    } else if (isTRUE(arg == 'Biweight')){

      for (i in 1:N){
        if (isTRUE(abs(x[i]) <= 1)){
          I <- 1
        }else {
          I <- 0
        }
        k[i] <- (15/16)* (1 - (x[i])^2)^2 * I
      }

    } else {

      for (i in 1:N){
        k[i] <- (1/sqrt(2 * pi))*exp(-(x[i]^2) / 2)
      }

    }

  } else {

#Scalar argument

    if (isTRUE(arg == 'Epanechnikov')){

      if (isTRUE(abs(x) <= 1)){
        I <- 1
      }else {
        I <- 0
      }

      k <- (3/4)*(1 - x^2)*I

    } else if (isTRUE(arg == 'Silverman')){

      k <- (1/2)*exp(-(abs(x)/sqrt(2)))*sin((abs(x)/sqrt(2))+(pi/4))

    } else if (isTRUE(arg == 'Biweight')){

      if (isTRUE(abs(x) <= 1)){
        I <- 1
      }else {
        I <- 0
      }

      k <- (15/16)* (1 - x^2)^2 * I

    } else {

      k <- (1/sqrt(2 * pi))*exp(-(x^2) / 2)

    }
  }
  k <- prod(k)
  return(k)
}
