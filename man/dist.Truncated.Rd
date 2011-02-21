\name{dist.Truncated}
\alias{dtrunc}
\alias{extrunc}
\alias{ptrunc}
\alias{qtrunc}
\alias{rtrunc}
\alias{vartrunc}
\title{Truncated Distributions}
\description{
  Density, distribution function, quantile function and random
  generation for truncated distributions.
}
\usage{
dtrunc(x, spec, a = -Inf, b = Inf, ...)
extrunc(spec, a = -Inf, b = Inf, ...)
ptrunc(x, spec, a = -Inf, b = Inf, ...)
qtrunc(p, spec, a = -Inf, b = Inf, ...)
rtrunc(n, spec, a = -Inf, b = Inf, ...)
vartrunc(spec, a = -Inf, b = Inf, ...)
}
\arguments{
  \item{n}{This is a the number of random draws for \code{rtrunc}.}
  \item{p}{This is a vector of probabilities.}
  \item{x}{This is a vector to be evaluated.}
  \item{spec}{The base name of a probability distribution is
    specified here. For example, to estimate the density of a
    truncated normal distribution, enter \code{norm}.}
  \item{a}{This is the lower bound of truncation, which defaults
    to negative infinity.}
  \item{b}{This is the upper bound of truncation, which defaults
    to infinity.}
  \item{\dots}{Additional arguments pertain to the probability
    distribution specified in the \code{spec} argument.}
}
\value{
  \code{dtrunc} gives the density, \code{extrunc} gives the expectation,
  \code{ptrunc} gives the distribution function, \code{qtrunc} gives the
  quantile function, \code{rtrunc} takes random draws, and
  \code{vartrunc} gives the variance of the truncated distribution.
}
\details{
  A truncated distribution is a conditional distribution that results
  from a priori restricting the domain of some other probability
  distribution. More than merely preventing values outside of truncated
  bounds, a proper truncated distribution integrates to one within the
  truncated bounds. In contrast to a truncated distribution, a censored
  distribution occurs when the probability distribution is still allowed
  outside of a pre-specified range.
}
\references{
Nadarajah, S. and Kotz, S. (2006). R Programs for Computing Truncated
Distributions. Journal of Statistical Software, 16, Code Snippet 2, p. 1--8.
}
\examples{
x <- seq(-3, 3, by = 0.1)
y <- dtrunc(x, "norm", a=-0.5, b=0.5, mean=0, sd=2)
}
\keyword{distribution}
