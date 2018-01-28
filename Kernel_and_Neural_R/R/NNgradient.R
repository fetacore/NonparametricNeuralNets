#'	The gradient of the three activation functions used in neural networks (sigmoidal,tanh and linear)
#'
#'	@param x	The inputs.
#'	@param arg The option of one of the three training functions (tanh by default)
#'	@return The	gradient of the activation function.
#'	@export

gradient <- function(x, arg='tanh'){

  if (isTRUE(arg == 'tanh')){
    g <- activation(x, 'tanh')
    z <- 1- g^2
  } else if (isTRUE(arg =='Sigmoid')){
    g <- activation(x, 'Sigmoid')
    z <- g * (1- g)
  } else {
    z <- matrix(data=1, nrow=NROW(x), ncol=NCOL(x))
  }
  return(z)
}
