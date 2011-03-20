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
    \item{Sigma}{This is a covariance matrix, and may be entered as
      either a matrix or vector.}
}
\value{
     \code{Cov2Cor} returns a correlation matrix, either as a matrix or
     vector according to whether the covariance matrix was passed to it
     as a matrix or vector.
}
\examples{
# The following example fits a model as Fit1, then adds a predictor, and
# fits another model, Fit2. The two models are compared with Bayes
# factors.

library(LaplacesDemon)
Sigma <- matrix(c(2,0.1,0.1,3),2,2)
Cov2Cor(Sigma)
Sigma.vector <- as.vector(Sigma)
Cov2Cor(Sigma.vector)
}
\keyword{covariance, correlation}