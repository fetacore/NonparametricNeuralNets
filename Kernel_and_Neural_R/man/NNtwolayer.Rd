% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNtwolayer.R
\name{NNtwolayer}
\alias{NNtwolayer}
\title{A three layer neural network for regression with quadratic loss.}
\usage{
NNtwolayer(y, x, iter = 100, hidden_1 = 20, hidden_2 = 50, eta = 0.005,
  ahidden_1 = "tanh", ahidden_2 = "Sigmoid", aout = "Linear")
}
\arguments{
\item{x}{The inputs (either vector or matrix).}

\item{iter}{The number of desired iterations.
@param hidden_1 The number of hidden units in the first hidden layer.
@param hidden_2 The number of hidden units in the second hidden layer.
@param eta The speed of learning (0.005 by default).
@param ahidden_1 The option of the activation function for the first hidden layer (tanh by default).
@param ahidden_2 The option of the activation function for the second hidden layer (sigmoid by default).}

\item{aout}{The option of the activation function for the output layer (linear by default).}
}
\value{
The output after training and the MSE per iteration (as a list).
}
\description{
@param y The target output.
}

