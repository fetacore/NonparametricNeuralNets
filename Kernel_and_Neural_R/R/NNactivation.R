#'	The three functions used in neural networks (sigmoidal, linear or tanh)
#'
#'	@param x	The inputs of the hidden neuron.
#'	@param arg The option of one of the two training functions (tanh by default)
#'	@return The	hidden neuron output z.
#'	@export


activation <- function(x, arg ='tanh'){
  if (isTRUE(arg == 'tanh')){
      z <- (exp(x)-exp(-x))/(exp(x)+exp(-x))
  } else if (isTRUE(arg =='Sigmoid')){
      z <- 1/(1+exp(-x))
  } else {
      z <- x
  }
  return(z)
}
