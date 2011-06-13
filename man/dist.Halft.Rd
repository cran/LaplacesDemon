\name{dist.Halft}
\alias{dhalft}
\alias{phalft}
\alias{qhalft}
\alias{rhalft}
\title{Half-t Distribution}
\description{
  These functions provide the density, distribution function, quantile function, and random generation for the half-t distribution.
}
\usage{
dhalft(x, scale=25, nu=1, log=FALSE)
phalft(q, scale=25, nu=1)
qhalft(p, scale=25, nu=1)
rhalft(n, scale=25, nu=1)
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{nu}{This is the scalar degrees of freedom parameter.}
  \item{scale}{This is the scale parameter alpha, which must be positive.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: p(theta) = (1 + (1/nu)*(theta/alpha)^2)^(-(nu+1)/2),
  theta >= 0
  \item Inventor: Derived from the Student t
  \item Notation 1: theta ~ Half-t(alpha, nu)
  \item Notation 2: theta ~ HT(alpha, nu)
  \item Notation 3: p(theta) = Half-t(theta | alpha, nu)
  \item Notation 4: p(theta) = HT(theta | alpha, nu)
  \item Parameter 1: scale parameter alpha > 0
  \item Parameter 2: degrees of freedom parameter nu
  \item Mean: E(theta) = unknown
  \item Variance: var(theta) = unknown
  \item Mode: 0
}

The Half-t distribution is an uninformative prior distribution for a
scale parameter that is more adaptable than the default recommended
half-Cauchy. When nu=1, the density is proportional to a proper half-Cauchy
distribution. When nu=-1, the density becomes an improper, uniform
prior distribution.
}
\value{
  \code{dhalft} gives the density,
  \code{phalft} gives the distribution function,
  \code{qhalft} gives the quantile function, and
  \code{rhalft} generates random deviates.
}
\examples{
x <- dhalft(1,25,1)
x <- phalft(1,25,1)
x <- qhalft(0.5,25,1)
x <- rhalft(10,25,1)
}
\keyword{distribution}