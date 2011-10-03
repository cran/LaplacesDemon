\name{dist.Laplace}
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
  \item{location}{This is the location parameter \eqn{\mu}{mu}.}
  \item{scale}{This is the scale parameter \eqn{\lambda}{lambda}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{1}{2 \lambda} \exp(-\frac{|\theta - \mu|}{\lambda})}{p(theta) = (1 / (2*lambda)) * exp(-(abs(theta - mu) / lambda))}
  \item Inventor: Pierre-Simon Laplace
  \item Notation 1: \eqn{\theta \sim}{theta ~} Laplace\eqn{(\mu,\lambda)}{(mu, lambda)}
  \item Notation 2: \eqn{\theta \sim \mathcal{L}(\mu, \lambda)}{theta ~ L(mu, lambda)}
  \item Notation 3: \eqn{p(\theta) =}{p(theta) =} Laplace\eqn{(\theta |
    \mu, \lambda)}{(theta | mu, lambda)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{L}(\theta | \mu, \lambda)}{p(theta) = L(theta | mu,lambda)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\lambda > 0}{lambda > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = 2 \lambda^2}{var(theta) = 2*lambda^2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The Laplace distribution is also called the double exponential 
distribution, because it looks like two exponential distributions back to
back with respect to location \eqn{\mu}{mu}. It is also called the
``First Law of Laplace'', just as the normal distribution is referred to
as the ``Second Law of Laplace''. The Laplace distribution is
symmetric with respect to \eqn{\mu}{mu}, though there are asymmetric
versions of the Laplace distribution. The PDF of the Laplace
distribution is reminiscent of the normal distribution; however,
whereas the normal distribution is expressed in terms of the squared
difference from the mean \eqn{\mu}{mu}, the Laplace density is
expressed in terms of the absolute difference from the mean,
\eqn{\mu}{mu}. Consequently, the Laplace distribution has fatter
tails than the normal distribution. It has been argued that the Laplace
distribution fits most things in nature better than the normal
distribution. Elsewhere, there are a large number of extensions to the
Laplace distribution, including asymmetric versions and
multivariate versions, among many more.

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