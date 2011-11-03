\name{MCSE}
\alias{MCSE}
\title{Monte Carlo Standard Error}
\description{
  Monte Carlo Standard Error (MCSE) is an estimate of the inaccuracy of
  Monte Carlo samples, usually regarding posterior samples from Markov
  chain Monte Carlo (MCMC) algorithms, such as with the
  \code{\link{LaplacesDemon}} function. MCSE approaches zero as the
  number of independent posterior samples approaches infinity. MCSE is
  essentially a standard deviation around the posterior mean of the
  samples, due to uncertainty associated with using an MCMC
  algorithm. The acceptable size of the MCSE depends on the acceptable
  uncertainty associated around the marginal posterior mean. Laplace's
  Demon prefers to continue updating until each MCSE is less than 6.7\%
  of each marginal posterior standard deviation (for more information on
  this stopping rule, see the \code{\link{Consort}} function).
}
\usage{
MCSE(x, method="IMPS", batch.size="sqrt", warn=TRUE)
}
\arguments{
  \item{x}{This is a vector of posterior samples for which MCSE will be
    estimated.}
  \item{method}{This is an optional argument for the method of
    estimation, and defaults to Geyer's \code{"IMPS"} method. Optional
    methods include \code{"sample.variance"} and \code{"batch.mean"}.
    Note that \code{"batch.mean"} is recommended only when the number of
    posterior samples is at least 1,000.}
  \item{batch.size}{This is an optional argument that corresponds only
    with \code{method="batch.means"}, and determines either the size of
    the batches (accepting a numerical argument) or the method of
    creating the size of batches, which is either \code{"sqrt"} or
    \code{"cuberoot"}, and refers to the length of \code{x}. The default
    argument is \code{"sqrt"}.}
  \item{warn}{Logical. If \code{warn=TRUE}, then a warning is provided
    with \code{method="batch.means"} whenever posterior sample size is
    less than 1,000.}
}
\details{
  The default method for estimating MCSE is Geyer's Initial Monotone
  Positive Sequence (IMPS) estimator (Geyer, 1992), which takes the
  asymptotic variance into account.

  The simplest method for estimating MCSE is to modify the formula for
  standard error, \eqn{\sigma(\textbf{x}) / \sqrt{N}}{sigma(x) /
    sqrt(N)}, to account for non-independence in the sequence
  \eqn{\textbf{x}}{x} of posterior samples. Non-independence is
  estimated with the \code{ESS} function for Effective Sample Size (see
  the \code{\link{ESS}} function for more details), where \eqn{M =
    ESS(\textbf{x})}{m = ESS(x)}, and MCSE is
  \eqn{\sigma(\textbf{x}) / \sqrt{M}}{sigma(x) / sqrt(M)}. Although this
  is the fastest and easiest method of estimation, it does not
  incorporate an estimate of the asymptotic variance of
  \eqn{\textbf{x}}{x}.

  The batch means method (Jones et al., 2006; Flegal et al., 2008)
  separates elements of \eqn{\textbf{x}}{x} into batches and estimates
  MCSE as a function of multiple batches. This method is excellent, but
  is not recommended when the number of posterior samples is less than
  1,000.
}
\references{
  Flegal, J.M., Haran, M., and Jones, G.L. (2008). "Markov chain Monte
  Carlo: Can We Trust the Third Significant Figure?". Statistical
  Science, 23, p. 250--260.

  Geyer, C.J. (1992). "Practical Markov Chain Monte Carlo". Statistical
  Science, 7, 4, p. 473--483.
  
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
  
  Jones, G.L., Haran, M., Caffo, B.S., and Neath, R. (2006). "Fixed-Width
  Output Analysis for Markov chain Monte Carlo". Journal of the American
  Statistical Association, 101, 1, p. 1537--1547.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{Consort}},
  \code{\link{ESS}},
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
MCSE(x)
MCSE(x, method="batch.means")
MCSE(x, method="sample.variance")
}
\keyword{mcse, utility}
