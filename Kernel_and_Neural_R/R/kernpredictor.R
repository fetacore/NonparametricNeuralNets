#' Non parametric predictor
#'
#' @param my_x The scalar value wrt which we do prediction
#' @param y   The given values before the prediction
#' @param x   The previous observations of x
#' @param h   The bandwidth
#' @param p     The order of the local polynomial (by default p=0)
#' @param distance The distance metric
#' @param kernel   The alphabetical argument for the kernel function of choice (from choose_kernel)
#' @return The predicted y corresponding to my_x
#' @export



kernpredictor <- function(my_x, y, x, h=0.5, p=0, distance='Eucledean', kernel='Gaussian'){
  require(MASS)
  T <- NROW(x)
  N <- NROW(my_x)
  W <- matrix(data = 0, nrow = T, ncol = T)

  k <- matrix(data=NA, nrow=1, ncol=T)
  for (j in 1:T){
    k[1,j] <- choice_of_kernel((distance_metric(my_x,x[j],distance)/h),kernel)
  }
  diag(W) <- k[1,]

  order = p + 1
  e <- matrix(data=1, nrow=order, ncol=1)
  B <- matrix(data=NA, nrow=T, ncol=order)
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
        B[j, o] <- (distance_metric(x[j],my_x,distance))^(o - 1)
      }
    }
  }
  f <- t(e) %*% ginv(t(B) %*% W %*% B) %*% t(B) %*% W %*% y
  return(f)
}
