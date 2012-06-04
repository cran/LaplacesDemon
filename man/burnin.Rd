\name{burnin}
\alias{burnin}
\title{Burn-in}
\description{
  The \code{burnin} function estimates the duration of burn-in in
  iterations for one or more Markov chains. ``Burn-in'' refers to the
  initial portion of a Markov chain that is not stationary and is still
  affected by its initial value.
}
\usage{
burnin(x, method="Geweke")
}
\arguments{
  \item{x}{This is a vector or matrix of posterior samples for which a
    the number of burn-in iterations will be estimated.}
  \item{method}{This argument defaults to \code{"Geweke"}, in which case
    stationarity is estimated with the \code{\link{Geweke.Diagnostic}}
    function. Alternatively, the \code{\link{KS.Diagnostic}} function
    may be used when \code{method="KS"} to estimate stationarity with
    the Kolmogorov-Smirnov test.}
}
\details{
  Burn-in is a colloquial term for the initial iterations in a Markov
  chain prior to its convergence to the target distribution. During
  burn-in, the chain is not considered to have ``forgotten'' its initial
  value.

  Burn-in is not a theoretical part of MCMC, but its use is the norm
  because of the need to limit the number of posterior samples due to
  computer memory. If burn-in were retained rather than discarded, then
  more posterior samples would have to be retained. If a Markov chain
  starts anywhere close to the center of its target distribution, then
  burn-in iterations do not need to be discarded.

  In the \code{\link{LaplacesDemon}} function, stationarity is estimated
  with the \code{\link{Geweke.Diagnostic}} function on all thinned
  posterior samples of each chain, beginning at cumulative 10\% intervals
  relative to the total number of samples, and the lowest number in
  which all chains are stationary is considered the burn-in.

  The term, ``burn-in'', originated in electronics regarding the initial
  testing of component failure at the factory to eliminate initial
  failures (Geyer, 2011).
}
\value{
  The \code{burnin} fucntion returns a vector equal in length to the
  number of MCMC chains in \code{x}, and each element indicates the
  maximum iteration in burn-in.
}
\references{
  Geyer, C.J. (2011). "Introduction to Markov Chain Monte Carlo". In
  S Brooks, A Gelman, G Jones, and M Xiao-Li (eds.), "Handbook of
  Markov Chain Monte Carlo", p. 3--48. Chapman \& Hall, Boca Raton, FL.

  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{Geweke.Diagnostic}},
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
burnin(x)
}
\keyword{diagnostic, stationarity, utility}
