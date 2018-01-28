#'	A three layer neural network for regression with quadratic loss.
#'
#'
#'  @param y The target output.
#'	@param x The inputs (either vector or matrix).
#'	@param iter The number of desired iterations.
#'  @param hidden_1 The number of hidden units in the first hidden layer.
#'  @param hidden_2 The number of hidden units in the second hidden layer.
#'  @param eta The speed of learning (0.005 by default).
#'  @param ahidden_1 The option of the activation function for the first hidden layer (tanh by default).
#'  @param ahidden_2 The option of the activation function for the second hidden layer (sigmoid by default).
#'	@param aout The option of the activation function for the output layer (linear by default).
#'	@return The output after training and the MSE per iteration (as a list).
#'	@importFrom stats runif
#'	@export

NNtwolayer <- function(y, x, iter=100, hidden_1=20, hidden_2=50, eta=0.005, ahidden_1='tanh', ahidden_2='Sigmoid', aout='Linear'){

  T <- NROW(x)
  R <- NCOL(x)
  inputs <- R
  outputs <- 1
  # Initial weights and biases with dimensions
  W_hidden_1 <- matrix(runif(inputs*hidden_1, 0, 1), nrow=inputs, ncol=hidden_1)
  B_hidden_1 <- matrix(1, nrow=1, ncol=hidden_1)
  W_hidden_2 <- matrix(runif(hidden_1*hidden_2, 0, 1), nrow=hidden_1, ncol=hidden_2)
  B_hidden_2 <- matrix(1, nrow=1, ncol=hidden_2)
  W_output <- matrix(runif(hidden_2, 0, 1), nrow=hidden_2, ncol=outputs)
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

      z_hidden_1 <- input %*% W_hidden_1 + B_hidden_1
      y_hidden_1 <- activation(z_hidden_1, ahidden_1)
      z_hidden_2 <- y_hidden_1 %*% W_hidden_2 + B_hidden_2
      y_hidden_2 <- activation(z_hidden_2, ahidden_2)
      z_output <- y_hidden_2 %*% W_output + B_output
      y_output[i] <- activation(z_output, aout)

      #Gradient descent

      d_output <- gradient(z_output, aout) * (y_output[i] - target)
      d_hidden_2 <-  (W_output %*% d_output) * t(gradient(z_hidden_2,ahidden_2))
      d_hidden_1 <-  (W_hidden_2 %*% d_hidden_2) * t(gradient(z_hidden_1,ahidden_1))

      gW_hidden_1 <- d_hidden_1 %*% input
      gW_hidden_2 <- d_hidden_2 %*% y_hidden_1
      gW_output <- d_output %*% y_hidden_2
      gB_hidden_1 <- d_hidden_1 * 1
      gB_hidden_2 <- d_hidden_2 * 1
      gB_output <- d_output * 1

      W_hidden_1 <- W_hidden_1 - eta * t(gW_hidden_1)
      w_hidden_2 <- W_hidden_2 - eta * t(gW_hidden_2)
      W_output <- W_output - eta * t(gW_output)
      B_hidden_1 <- B_hidden_1 - eta * t(gB_hidden_1)
      B_hidden_2 <- B_hidden_2 - eta * t(gB_hidden_2)
      B_output <- B_output - eta * t(gB_output)

      #Error

      error[i] <- 0.5 * (y_output[i] - target)^2
    }
    print(j)
    mse[j] <- mean(error)
    j <- j + 1
  }
#  return(list(y_output , mse))
  return(y_output)
}
