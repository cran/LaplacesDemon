\name{dist.Categorical}
\alias{dcat}
\alias{rcat}
\title{The Categorical Distribution}
\description{
  This is the density and random deviates function for the categorical
  distribution with probabilities parameter \eqn{p}.
}
\usage{
dcat(x, p, log=FALSE)
rcat(n, p)
}
\arguments{
  \item{x}{This is a vector of discrete data with \eqn{k} discrete
    categories, and is of length \eqn{n}. This function also accepts
    \eqn{x} after it has been converted to an \eqn{n \times k}{n x k}
    indicator matrix, such as with the \code{as.indicator.matrix} function.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{p}{This is a vector of length \eqn{k} or \eqn{n \times k}{n x k}
    matrix of probabilities.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
    \itemize{
  \item Application: Discrete Univariate
  \item Density: \eqn{p(\theta) = \sum \theta p}{p(theta) = Sum (theta * p)}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{CAT}(p)}{theta ~ CAT(p)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{CAT}(\theta | p)}{p(theta) = CAT(theta | p)}
  \item Parameter 1: probabilities \eqn{p}
  \item Mean: \eqn{E(\theta)}{E(theta)} = Unknown
  \item Variance: \eqn{var(\theta)}{var(theta)} = Unknown
  \item Mode: \eqn{mode(\theta)}{mode(theta)} = Unknown
}
Also called the discrete distribution, the categorical distribution
describes the result of a random event that can take on one of \eqn{k}
possible outcomes, with the probability \eqn{p} of each outcome
separately specified. The vector \eqn{p} of probabilities for each
event must sum to 1. The categorical distribution is often used, for
example, in the multinomial logit model. The conjugate prior is the
Dirichlet distribution.
}
\value{
  \code{dcat} gives the density and
  \code{rcat} generates random deviates.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{as.indicator.matrix}},
  \code{\link{ddirichlet}},
  \code{\link{dmultinom}}.
}
\examples{
library(LaplacesDemon)
dcat(x=1, p=c(0.3,0.3,0.4))
rcat(n=10, p=c(0.1,0.3,0.6))
}
\keyword{distribution}
