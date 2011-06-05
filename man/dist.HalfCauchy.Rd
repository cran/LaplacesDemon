\name{dist.HalfCauchy}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{qhalfcauchy}
\alias{rhalfcauchy}
\title{Half-Cauchy Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the half-Cauchy distribution.
}
\usage{
dhalfcauchy(x, scale=25, log=FALSE)
phalfcauchy(q, scale=25)
qhalfcauchy(p, scale=25)
rhalfcauchy(n, scale=25)
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{scale}{This is the scale parameter alpha, which must be positive.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: p(theta) = 2alpha / pi(theta^2 + alpha^2), theta >= 0
  \item Inventor: Derived from Cauchy
  \item Notation 1: theta ~ Half-Cauchy(alpha)
  \item Notation 2: theta ~ HC(alpha)
  \item Notation 3: p(theta) = Half-Cauchy(theta | alpha)
  \item Notation 4: p(theta) = HC(theta | alpha)
  \item Parameter 1: scale parameter alpha > 0
  \item Mean: E(theta) = does not exist
  \item Variance: var(theta) = does not exist
  \item Mode: 0
}

The Half-Cauchy distribution with scale alpha=25 is a recommended
default uninformative prior distribution for a standard deviation parameter.
}
\value{
  \code{dhalfcauchy} gives the density,
  \code{phalfcauchy} gives the distribution function,
  \code{qhalfcauchy} gives the quantile function, and
  \code{rhalfcauchy} generates random deviates.
}
\examples{
x <- dhalfcauchy(1,25)
x <- phalfcauchy(1,25)
x <- qhalfcauchy(0.5,25)
x <- rhalfcauchy(1,25)
}
\keyword{distribution}