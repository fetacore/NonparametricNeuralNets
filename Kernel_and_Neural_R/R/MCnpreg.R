#' A Monte Carlo simulation
#'
#' @param true The true model
#' @param error The error component of the data generating function
#' @param x The regressor(s)
#' @param testh If we want to test the evolution of the CV criterion with different smoothing parameters (False by default)
#' @param p The order of the local polynomial fit (local constant by default)
#' @param method The kernel function (Gaussian by default)
#' @export

MCnpreg <- function(true, error, x, testh=FALSE, p=0, method = 'Gaussian'){

  y <- true + error
  T <- NROW(x)
  if (isTRUE(testh==FALSE)){
    Bias <- matrix(data=NA, nrow=T, ncol=1)
    Var <- matrix(data=NA, nrow=T, ncol=1)
    MISE <- 1
    h <- the_bandwidth(x)
    f <- nonparpoly(y,x,h,p,method)
    Bias <- true - f
    Var <- (sd(true - f))^2
    MISE <- mean((true - f)^2)
  } else {
    h <- runif(10000, 0.001, 50)
    A <- NROW(h)
    Bias <- matrix(data=NA, nrow=T, ncol=A)
    Var <- matrix(data=NA, nrow=T, ncol=A)
    ovBias <- matrix(data=NA, nrow=A, ncol=1)
    ovVar <- matrix(data=NA, nrow=A, ncol=1)
    MISE <- matrix(data=NA, nrow=A, ncol=1)
    for (a in 1:NROW(h)){
      f <- nonparpoly(y,x,h[a],p,method)
      Bias[,a] <- true - f
      Var[,a] <- (sd(true - f))^2
      ovBias[a] <- sum(true-f)
      ovVar[a] <- sum((sd(true - f))^2)
      MISE[a] <- mean((true - f)^2)
      print(a)
    }
    return(c(h, Bias, Var, MISE, ovBias, ovVar))
  }









}
