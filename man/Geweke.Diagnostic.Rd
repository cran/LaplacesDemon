\name{Geweke.Diagnostic}
\alias{Geweke.Diagnostic}
\title{Geweke's Convergence Diagnostic}
\description{
  Geweke (1992) proposed a convergence diagnostic for Markov chains.
  This diagnostic is based on a test for equality of the means of the
  first and last part of a Markov chain (by default the first 10\% and
  the last 50\%).  If the samples are drawn from a stationary
  distribution of the chain, then the two means are equal and Geweke's
  statistic has an asymptotically standard normal distribution.

  The test statistic is a standard Z-score: the difference between the
  two sample means divided by its estimated standard error.  The
  standard error is estimated from the spectral density at zero, and so
  takes into account any autocorrelation.

  The Z-score is calculated under the assumption that the two parts of
  the chain are asymptotically independent.
}
\usage{Geweke.Diagnostic(x)}
\arguments{
     \item{x}{
          This required argument is a vector or matrix of posterior
          samples.
          }
}
\details{
  The \code{Geweke.Diagnostic} is essentially the same as the
  \code{geweke.diag} function in the \code{coda} package, but
  programmed to accept a simple vector or matrix, so it does not require
  an \code{mcmc} object.
}
\value{
  A vector is returned, in which each element is a Z-score for a test of
  equality that compares the means of the first and last parts of each
  chain supplied as \code{x} to \code{Geweke.Diagnostic}.
}
\references{
  Geweke, J. (1992). Evaluating the Accuracy of Sampling-Based
  Approaches to Calculating Posterior Moments. In
  \emph{Bayesian Statistics 4} (ed JM Bernardo, JO Berger, AP Dawid,
  and AFM Smith). Clarendon Press, Oxford, UK.
}
\seealso{\code{\link{LaplacesDemon}}}
