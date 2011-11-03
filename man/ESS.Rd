\name{ESS}
\alias{ESS}
\title{Effective Sample Size due to Autocorrelation}
\description{
  This function may be used to estimate the effective sample size (ESS)
  of a target distribution, where the sample size is reduced by
  autocorrelation. ESS is a measure of how well the chains are mixing.
}
\usage{ESS(x)}
\arguments{
     \item{x}{
          This required argument is a vector or matrix of posterior
          samples.
          }
}
\details{
  Effective Sample Size (ESS) was recommended by Radford Neal in the
  panel discussion of Kass et al. (1998). The \code{ESS} function
  uses the \code{ar} function to estimate autocorrelation. When the
  posterior distribution is sampled with a Markov chain Monte Carlo
  (MCMC) algorithm, there is usually autocorrelation present in the
  samples. More autocorrelation is associated with less posterior
  sampled information, because the information in the samples is
  autocorrelated, or put another way, successive samples are not
  independent from earlier samples. This reduces the effective sample
  size of, and precision in representing, the posterior
  distribution. \code{ESS} is one of the criteria in the
  \code{\link{Consort}} function, where stopping the MCMC updates is
  not recommended until \code{ESS} \eqn{\ge 100}. Although the need
  for precision of each modeler differs with each model, it is often
  a good goal to obtain \code{ESS} \eqn{= 1000}.

  This is numerically the same as the \code{effectiveSize} function in
  the \code{coda} package, but programmed to accept a simple vector or
  matrix so it does not require an \code{mcmc} or \code{mcmc.list}
  object, and the result is bound to be less than or equal to the
  original number of samples.
}
\value{
  A vector is returned, and each element is the effective sample size
  (ESS) for a corresponding column of \code{x}, after autocorrelation has
  been taken into account.
}
\references{
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}

  Kass, R.E., Carlin, B.P., Gelman, A., and Neal, R. (1998). Markov
  Chain Monte Carlo in Practice: A Roundtable Discussion. The American
  Statistician, 52, p. 93--100. 
}
\seealso{
  \code{\link{Consort}},
  \code{\link{LaplacesDemon}}.
}
