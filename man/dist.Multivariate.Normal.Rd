\name{dist.Multivariate.Normal}
\alias{dmvn}
\alias{rmvn}
\title{Multivariate Normal Distribution}
\description{
  These functions provide the density and random number generation for
  the multivariate normal distribution.
}
\usage{
dmvn(x, mu, Sigma, log=FALSE) 
rmvn(n=1, mu, Sigma)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
    \code{k} or a matrix with \code{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is the mean vector with length \code{k} or matrix with
    \code{k} columns. The \code{rmvn} function requires a vector.}
  \item{Sigma}{This is the \code{k x k} covariance matrix.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: p(theta) = (1/((2*pi)^(k/2)*|Sigma|^(1/2)))*exp(-.5*(theta-mu)'*Sigma^(-1)*(theta-mu))
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: theta ~ MVN(mu, Sigma)
  \item Notation 2: theta ~ N[k](mu, Sigma)
  \item Notation 3: p(theta) = MVN(theta | mu, Sigma)
  \item Notation 4: p(theta) = N[k](theta | mu, Sigma)
  \item Parameter 1: location vector mu
  \item Parameter 2: positive-definite covariance matrix Sigma
  \item Mean: E(theta) = mu
  \item Variance: var(theta) = Sigma
  \item Mode: mu
}

The multivariate normal distribution, or multivariate Gaussian
distribution, is a multidimensional extension of the one-dimensional
or univariate normal or Gaussian distribution. A random vector is
considered to be multivariate normally distributed if every linear
combination of its components has a univariate normal distribution.
This distribution has a mean parameter vector \code{mu} of length
\code{k}, and a \code{k x k} covariance matrix \code{Sigma}, which
must be positive-definite.

These functions are similar to those in the \code{bayesm}, \code{mnormt}, and
\code{mvtnorm} packages.
}
\value{
  \code{dmvn} gives the density, and 
  \code{rmvn} generates random deviates.
}
\examples{
x <- dmvn(c(1,2,3), c(0,1,2), matrix(c(1,0,0,0,1,0,0,0,1),3,3))
x <- rmvn(100, c(0,1,2), matrix(c(1,0,0,0,1,0,0,0,1),3,3))
}
\keyword{distribution}