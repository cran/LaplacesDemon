\name{as.initial.values}
\alias{as.initial.values}
\title{Initial Values}
\description{
  This function returns the most recent posterior samples from an object
  of class \code{demonoid}, the most recent posterior samples from the
  chain with the lowest deviance in an object of class
  \code{demonoid.hpc}, the posterior modes of an object of class
  \code{laplace}, the posterior means of an object of class \code{pmc}
  with one mixture component, or the latest means of the importance
  sampling distribution of an object of class \code{pmc} with multiple
  mixture components. The returned values are intended to be the initial
  values for future updates.
}
\usage{
as.initial.values(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid},
    \code{demonoid.hpc}, \code{laplace}, or \code{pmc}.}
}
\details{
  Unless it is known beforehand how many iterations are required for
  Laplace Approximation to converge, MCMC to appear converged, or the
  normalized perplexity to stabilize in PMC, multiple updates are
  necessary. An additional update, however, should not begin with the
  same initial values as the original update, because it will have to
  repeat the work already accomplished. For this reason, the
  \code{as.initial.values} function may be used at the end of an update
  to change the previous initial values to the latest values.
}
\value{
  The returned value is a vector (or matrix in the case of PMC with
  multiple mixture components) of the latest values, which may now be
  used as initial values for a future update.
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}}, and
  \code{\link{PMC}}.
}
\keyword{Initial Values}