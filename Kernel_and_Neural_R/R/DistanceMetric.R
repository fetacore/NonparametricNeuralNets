#'	Defines the desired metric for our problem.
#'
#'	@param x	A	number/vector.
#'	@param y	A	number/vector.
#'	@param arg The optional distance metric (Alphabetical/Eucledean by default)
#'	@param p For the Lp-norm (requires arg='L')
#'	@note If argument is not alphabetical, it returns the Eucledean Metric
#'	@return The	distance	of	x	and	y.
#'	@export



distance_metric <- function(x, y, arg='Eucledean', p=NULL){
  o <- rep(0, NROW(x))

  if (isTRUE(arg=='Eucledean')){
    for (i in 1:NROW(x)){
      o[i] <- (x[i] - y[i])^2
    }
    omega <- sqrt(sum(o))

  } else {
    if (arg == 'L'){
      for (i in 1:NROW(x)){
        o[i] <- abs((x[i] - y[i]))^p
      }
      omega <- (sum(o))^(1/p)

    }else if (arg == 'exp'){
      for (i in 1:NROW(x)){
        o[i] <- abs(exp(x[i]) - exp(y[i]))
      }
      omega <- sum(o)

    } else if (arg == 'max'){
      for (i in 1:NROW(x)){
        o[i] <- abs(x[i] - y[i])
      }
      t <- which.max(o)
      omega <- o[t]

    } else if (arg == 'Mahalanobis'){
      S <- cov(x,y)
      omega <- sqrt(t(x-y) %*% solve(S) %*% (x-y))

    } else {
      for (i in 1:NROW(x)){
        o[i] <- abs(x[i] - y[i])
      }
      omega <- sum(o)
    }
  }
  return(omega)
}
