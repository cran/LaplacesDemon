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
  \item{p}{This is a vector of probabilities.}
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
  \item Density: \eqn{p(\theta) = (1 + \frac{1}{\nu} (\theta /
    \alpha)^2)^{(-(\nu+1)/2)}, \quad \theta \ge 0}{p(theta) = (1 + (1/nu)*(theta/alpha)^2)^(-(nu+1)/2), theta >= 0}
  \item Inventor: Derived from the Student t
  \item Notation 1: \eqn{\theta \sim \mathcal{HT}(\alpha, \nu)}{theta ~ HT(alpha, nu)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{HT}(\theta | \alpha, \nu)}{p(theta) = HT(theta | alpha, nu)}
  \item Parameter 1: scale parameter \eqn{\alpha > 0}{alpha > 0}
  \item Parameter 2: degrees of freedom parameter \eqn{\nu}{nu}
  \item Mean: \eqn{E(\theta)}{E(theta)} = unknown
  \item Variance: \eqn{var(\theta)}{var(theta)} = unknown
  \item Mode: \eqn{mode(\theta) = 0}{mode(theta) = 0}
}

The half-t distribution is an uninformative prior distribution for a
scale parameter that is more adaptable than the default recommended
half-Cauchy. When \eqn{\nu=1}{nu=1}, the density is proportional to
a proper half-Cauchy distribution. When \eqn{\nu=-1}{nu=-1}, the
density becomes an improper, uniform prior distribution.
}
\value{
  \code{dhalft} gives the density,
  \code{phalft} gives the distribution function,
  \code{qhalft} gives the quantile function, and
  \code{rhalft} generates random deviates.
}
\seealso{
  \code{\link{dst}} for the Student-t distribution.
}
\examples{
x <- dhalft(1,25,1)
x <- phalft(1,25,1)
x <- qhalft(0.5,25,1)
x <- rhalft(10,25,1)
}
\keyword{distribution}