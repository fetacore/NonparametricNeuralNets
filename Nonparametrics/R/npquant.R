#' Nonparametric Quantile Regression
#'
#' @param y The regressant
#' @param x The regressor
#' @param a Which quantile
#' @param hx The bandwidth for xs
#' @param hy The bandwidth for ys
#' @param method For the kernel function
#' @export

npquant <- function(y,x,a=0.5,hx=0.05,hy=0.1,method='Gaussian'){

  acumulated.distrib= function(sample,x){
    minors= 0
    for(n in sample){
      if(n<=x){
        minors= minors+1
      }
    }
    return (minors/length(sample))
  }
  T <- NROW(x)
  k <- kernmatrix(x,hx,method)
  f <- parzen_rosenblatt(x,hx,method)
  Fhat <- matrix(data=0,nrow=T,ncol=T)
  G <- matrix(data=0,nrow=T,ncol=T)
  q <- matrix(data=0,nrow=T,ncol=1)
  Fhatinv <- matrix(data=0,nrow=T,ncol=T)
  for (i in 1:T){
    for (e in 1:T){
      yp <- (y[i]-y[e])/hy
      G[i,e] <- acumulated.distrib(y,yp)
    }
    for (j in 1:T){
      Fhat[i,j] <- ((1/T)*sum(k[j,]*G[i,]))/(f[i])
      if (isTRUE(Fhat[i,j] >a) | isTRUE(Fhat[i,j]==a)){
        Fhatinv[i,j] <- Fhat[i,j]
      } else {
        Fhatinv[i,j] <- 1
      }
    }
    q[i] <- y[which.min(Fhatinv[i,])]
  }
  return(list(q,G,Fhat,Fhatinv))
}
