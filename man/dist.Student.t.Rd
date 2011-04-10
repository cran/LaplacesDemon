\name{dist.Student.t}
\alias{dst}
\alias{pst}
\alias{qst}
\alias{rst}
\title{Student t Distribution: Univariate}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate Student t
  distribution with location parameter \code{mu}, scale parameter
  \code{sigma}, and degrees of freedom parameter \code{nu}.
}
\usage{
dst(x, mu=0, sigma=1, nu=10, log=FALSE)
pst(q, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
qst(p, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
rst(n, mu=0, sigma=1, nu=10)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter mu.}
  \item{sigma}{This is the scale parameter sigma, which must be
    positive.}
  \item{nu}{This is the degrees of freedom parameter nu, which must be
    positive.}
  \item{lower.tail}{Logical. If TRUE, then probabilities are P[X <= x],
          otherwise, P[X > x].}
  \item{log, log.p}{Logical. If \code{log=TRUE} then the logarithm of the
    density or probability is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: p(theta) = [{Gamma[(nu + 1) / 2]} / Gamma(nu/2)]*sqrt(nu*pi)*sigma*[1 + (1/nu)*[(theta - mu)/sigma]^2]^[(-nu + 1)/2]
  \item Inventor: William Sealy Gosset (1908)
  \item Notation 1: theta ~ t(mu, sigma, nu)
  \item Notation 2: p(theta) = t(theta | mu, sigma, nu)
  \item Parameter 1: location parameter mu
  \item Parameter 2: scale parameter sigma > 0
  \item Parameter 3: degrees of freedom nu > 0
  \item Mean: E(theta) = mu, for nu > 1, otherwise undefined
  \item Variance: var(theta) = [nu / (nu - 2)]sigma^2, for nu > 2
  \item Mode: mu
}

The Student t-distribution is often used as an alternative to the normal distribution as a model for data. It is frequently the case that real data have heavier tails than the normal distribution allows for. The classical approach was to identify outliers and exclude or downweight them in some way. However, it is not always easy to identify outliers (especially in high dimensions), and the Student t-distribution is a natural choice of model-form for such data and provides a parametric approach to robust statistics.

The degrees of freedom parameter, nu, controls the kurtosis of the
distribution and is correlated with the scale parameter sigma. The
likelihood can have multiple local maxima and, as such, it is often
necessary to fix nu at a fairly low value and estimate the other
parameters taking this as given. Some authors report that values between
3 and 9 are often good choices, and some authors suggest 5 is often a
good choice.

In the limit nu -> infinity, the Student t-distribution approaches N(mu, sigma^2). The case of nu = 1 is the Cauchy distribution.

These functions are similar to those in the \code{gamlss.dist} package.
}
\value{
  \code{dst} gives the density,
  \code{pst} gives the distribution function,
  \code{qst} gives the quantile function, and
  \code{rst} generates random deviates.
}
\examples{
x <- dst(1,0,1,10)
x <- pst(1,0,1,10)
x <- qst(0.5,0,1,10)
x <- rst(100,0,1,10)
}
\keyword{distribution}