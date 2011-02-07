\name{laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
\title{Laplace Distribution: Univariate Symmetric}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, symmetric, Laplace
  distribution with location parameter \code{location} and scale
  parameter \code{scale}.
}
\usage{
dlaplace(x, location=0, scale=1, log=FALSE)
plaplace(q, location=0, scale=1)
qlaplace(p, location=0, scale=1)
rlaplace(n, location=0, scale=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{location}{This is the location parameter mu.}
  \item{scale}{This is the scale parameter b, which must be positive.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: p(theta) = (1 / (2*b)) * exp(-(abs(theta - mu) / b))
  \item Inventor: Pierre-Simon Laplace
  \item Notation 1: theta ~ Laplace(mu, b)
  \item Notation 2: theta ~ L(mu, b)
  \item Notation 3: p(theta) = Laplace(theta | mu, b)
  \item Notation 4: p(theta) = L(theta | mu,b)
  \item Parameter 1: location parameter mu
  \item Parameter 2: scale parameter b > 0
  \item Mean: E(theta) = mu
  \item Variance: var(theta) = 2*b^2
  \item Mode: mu
}

The Laplace distribution is also called the double exponential 
distribution, because it looks like two exponential distributions back to
back with respect to location mu. The PDF of the Laplace distribution is
reminiscent of the normal distribution; however, whereas the normal
distribution is expressed in terms of the squared difference from the mean
mu, the Laplace density is expressed in terms of the absolute difference
from the mean, mu. Consequently, the Laplace distribution has fatter
tails than the normal distribution. It has been argued that the Laplace
distribution fits most things in nature better than the normal
distribution. Elsewhere, there are also asymmetric versions and
multivariate versions of the Laplace distribution.

These functions are similar to those in the \code{VGAM} package.
}
\value{
  \code{dlaplace} gives the density,
  \code{plaplace} gives the distribution function,
  \code{qlaplace} gives the quantile function, and
  \code{rlaplace} generates random deviates.
}
\examples{
x <- dlaplace(1,0,1)
x <- plaplace(1,0,1)
x <- qlaplace(0.5,0,1)
x <- rlaplace(100,0,1)
}
\keyword{distribution}