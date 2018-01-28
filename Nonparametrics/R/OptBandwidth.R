#' Automated Bandwidth Selection
#'
#' @param x   The observations
#' @param method    The desired method of selection (Silverman's Rule of Thumb by default, if numerical then this is the user defined value of h)
#' @param y   If we need the cross validation criterion the regressant
#' @param p   For the nonparpoly function (order of polynomial local fit)
#' @param k   If we need the k-nearest neighborhood
#' @return The resulting bandwidth (either scalar or vector)
#' @importFrom stats IQR sd optim
#' @export

the_bandwidth <- function(x, method='Silverman',y=NULL, p=1, k=NULL){
  T <- NROW(x)


  if(isTRUE(method == 'Silverman')){
    h <- (1/T)^(1/5) * sd(x) * 1.06
  } else {

    if(isTRUE(method == 'Robust')){

      h <- 1.06 * min( sd(x) , IQR(x)/1.34) * (1/T)^(1/5)

    }  else {
      h <- method
    }

  }

  return(h)
}
