\name{Cov2Cor}
\alias{Cov2Cor}
\title{Covariance to Correlation}
\description{
This function calculates a correlation matrix from a covariance matrix,
and accepts the covariance matrix either in matrix or vector form. This
function may be useful inside a model specification and also with
converting posterior draws of the elements of a covariance matrix to a
correlation matrix.
}
\usage{
Cov2Cor(Sigma)
}
\arguments{
     \item{Sigma}{This is a covariance matrix, \eqn{\Sigma}{Sigma}, and
          may be entered as either a matrix or vector.}
}
\details{\code{Cov2Cor} is an expanded form of the \code{cov2cor} function in the \code{stats} package, where \code{Cov2Cor} is also able to accept and return a vectorized matrix.}
\value{
     \code{Cov2Cor} returns a correlation matrix, either as a matrix or
     vector according to whether the covariance matrix was passed to it
     as a matrix or vector.
}
\examples{
library(LaplacesDemon)
Sigma <- matrix(c(2,0.1,0.1,3),2,2)
Cov2Cor(Sigma)
Sigma.vector <- as.vector(Sigma)
Cov2Cor(Sigma.vector)
}
\keyword{covariance, correlation}