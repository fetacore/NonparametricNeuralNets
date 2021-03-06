% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernmatrix.R
\name{kernmatrix}
\alias{kernmatrix}
\title{The Kernel Matrix}
\usage{
kernmatrix(x, h, kernel = "Gaussian")
}
\arguments{
\item{x}{The observations (vector or matrix)}

\item{h}{The bandwidth}

\item{kernel}{The alphabetical argument for the kernel function of choice (from choose_kernel)}
}
\value{
The k((x[i]-x[j])/h) matrix with i rows and j columns
}
\description{
The Kernel Matrix
}
\note{
By default we use the Gaussian Kernel
}

