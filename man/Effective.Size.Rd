\name{Effective.Size}
\alias{Effective.Size}
\title{Effective Sample Size due to Autocorrelation}
\description{
  This may be used to estimate the effective sample size of a target
  distribution, where the sample size is reduced by autocorrelation.
}
\usage{Effective.Size(x)}
\arguments{
     \item{x}{
          This required argument is a vector or matrix of posterior
          samples.
          }
}
\details{
The \code{Effective.Size} function uses the \code{ar} function to
estimate autocorrelation.

This is numerically the same as the \code{effectiveSize} function in
the \code{coda} package, but programmed to accept a simple vector or
matrix, so it does not require an \code{mcmc} or \code{mcmc.list}
object.
}
\value{
  A vector is returned, and each element is the effective sample size
  for a corresponding column of \code{x}, after autocorrelation has
  been taken into account.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
