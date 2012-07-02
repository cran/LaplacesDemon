\name{as.initial.values}
\alias{as.initial.values}
\title{Initial Values}
\description{
  This function returns the most recent posterior samples from an object
  of class \code{demonoid}, the most recent posterior samples from the
  chain with the lowest deviance in an object of class
  \code{demonoid.hpc}, or the posterior modes of an object of class
  \code{laplace}. The returned values are intended to be the initial
  values for future updates.
}
\usage{
as.initial.values(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid},
    \code{demonoid.hpc}, or \code{laplace}.}
}
\details{
  Unless it is known beforehand how many iterations are required for
  Laplace Approximation to converge or MCMC to appear converged,
  multiple updates are necessary. An additional update, however, should
  not begin with the same initial values as the original update, because
  it will have to repeat the work already accomplished. For this reason,
  the \code{as.initial.values} function may be used at the end of an
  update to change the previous initial values to the latest values. 
}
\value{
     The returned value is a vector of the latest values, which may now
     be used as initial values for a future update.
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{LaplacesDemon.hpc}}.
}
\keyword{Initial Values}