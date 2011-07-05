\name{dist.Categorical}
\alias{dcat}
\title{The Categorical Distribution}
\description{
  This is the density function for the categorical distribution with parameter \code{p}.
}
\usage{
dcat(x, p, log=FALSE)
}
\arguments{
  \item{x}{vector of discrete data with \code{k} discrete categories,
    and of length \code{n}. This function also accepts \code{x} after it
  has been converted to an \code{n x k} indicator matrix.}
  \item{p}{vector (of length \code{k}) or \code{n x k} matrix of probabilities.}
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{dcat} gives the density for the categorical distribution.
}
\details{
    \itemize{
  \item Application: Discrete Univariate
  \item Density: p(theta) = Sum (theta * p)
  \item Inventor: Unknown
  \item Notation 1: theta ~ Cat(p)
  \item Notation 2: p(theta) = Cat(theta | p)
  \item Parameter 1: probabilities p
  \item Mean: E(theta) = Unknown
  \item Variance: var(theta) = Unknown
  \item Mode: Unknown
}
Also called the discrete distribution, the categorical distribution
describes the result of a random event that can take on one of \code{k}
possible outcomes, with the probability \code{p} of each outcome
separately specified. The vector \code{p} of probabilities for each
event must sum to 1. The categorical distribution is often used, for
example, in the multinomial logit model. The conjugate prior is the
Dirichlet distribution.
}
\examples{
dcat(1, c(0.3,0.3,0.4))
}
\keyword{distribution}
