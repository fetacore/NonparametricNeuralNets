#'	A two layer neural network for regression with quadratic loss.
#'
#'
#'  @param y The target output.
#'	@param x The inputs (either vector or matrix).
#'	@param iter The number of desired iterations.
#'  @param hidden The number of hidden units in the hidden layer.
#'  @param eta The speed of learning (0.005 by default).
#'  @param ahidden The option of the activation function for the hidden layer (tanh by default).
#'	@param aout The option of the activation function for the output layer (linear by default).
#'	@return The output after training and the MSE per iteration (as a list).
#'	@importFrom stats runif
#'	@export

NNonelayer <- function(y, x, iter=100, hidden=20, eta=0.005, ahidden='tanh', aout='Linear'){

  T <- NROW(x)
  R <- NCOL(x)
  inputs <- R
  outputs <- 1
  # Initial weights and biases with dimensions
  W_hidden <- matrix(runif(inputs*hidden, 0, 1), nrow=inputs, ncol=hidden)
  B_hidden <- matrix(1, nrow=1, ncol=hidden)
  W_output <- matrix(runif(hidden, 0, 1), nrow=hidden, ncol=outputs)
  B_output <- matrix(1, nrow=1, ncol=outputs)
  y_output <- matrix(data=NA, nrow=T, ncol=1)
  mse <- matrix(data=0, nrow=1, ncol=iter)
  j <- 1
  while (!isTRUE(j == iter+1)){
    error <- matrix(data=0, nrow=T, ncol=1)
    for (i in 1:T){
      target <- y[i]

      if (isTRUE(R==1)){
        input <- (x[i]-mean(x))/sd(x)
      } else {
        input <- matrix(data=NA, nrow=1, ncol=R)
        for (r in 1:R){
          input[r] <- (x[i,r]-mean(x[,r]))/sd(x[,r])
        }
      }

      #The Backpropagation algorithm

      z_hidden <- input %*% W_hidden + B_hidden
      y_hidden <- activation(z_hidden, ahidden)
      z_output <- y_hidden %*% W_output + B_output
      y_output[i] <- activation(z_output, aout)

      #Gradient descent

      d_output <- gradient(z_output, aout) * (y_output[i] - target)
      d_hidden <- (t(gradient(z_hidden, ahidden)) %*% d_output) * W_output

      gW_hidden <- d_hidden %*% input
      gW_output <- t(d_output) %*% y_hidden
      gB_hidden <- d_hidden * 1
      gB_output <- d_output * 1

      W_hidden <- W_hidden - eta * t(gW_hidden)
      W_output <- W_output - eta * t(gW_output)
      B_hidden <- B_hidden - eta * t(gB_hidden)
      B_output <- B_output - eta * t(gB_output)

      #Error

      error[i] <- 0.5 * (y_output[i] - target)^2
    }
    print(j)
    mse[j] <- mean(error)
    j <- j + 1
  }
  return(list(y_output , mse))
}
