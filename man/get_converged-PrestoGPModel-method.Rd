% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrestoGP_Model.R
\name{get_converged,PrestoGPModel-method}
\alias{get_converged,PrestoGPModel-method}
\alias{get_converged}
\title{Did the PrestoGP model converge?}
\usage{
\S4method{get_converged}{PrestoGPModel}(model)
}
\arguments{
\item{model}{The PrestoGP model object}
}
\value{
A boolean value indicating whether or not the model converged.
}
\description{
This method returns a boolean value that indicates whether or not a
PrestoGP model converged.
}
\examples{
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
y <- soil[,4]                   # predict moisture content
X <- as.matrix(soil[,5:9])
locs <- as.matrix(soil[,1:2])

soil.vm <- new("VecchiaModel", n_neighbors = 10)
soil.vm <- prestogp_fit(soil.vm, y, X, locs)
get_converged(soil.vm)
}
\references{
\itemize{
\item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
land-use regression for ground-level nitrogen dioxide", The Annals of
Applied Statistics (2021) 15(2):688-710.
}
}
\seealso{
\code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
}
