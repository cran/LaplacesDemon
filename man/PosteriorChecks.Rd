\name{PosteriorChecks}
\alias{PosteriorChecks}
\title{Posterior Checks}
\description{
  Not to be confused with posterior predictive checks, this function
  provides additional information about the marginal posterior
  distributions, such as the probability that each posterior coefficient
  of the parameters (referred to generically as \eqn{\theta}{theta}), is
  greater than zero [\eqn{p(\theta > 0)}{p(theta > 0)}], the estimated
  number of modes, the kurtosis and skewness of the posterior
  distributions, the burn-in of each chain (for MCMC only), and the
  integrated autocorrelation time. A correlation matrix is
  provided only for objects of class \code{demonoid}.
}
\usage{
PosteriorChecks(x, Parms)
}
\arguments{
  \item{x}{This required argument accepts an object of class
    \code{demonoid} or \code{laplace}.}
  \item{Parms}{This argument accepts a vector of quoted strings to be
    matched for selecting parameters. This argument defaults to
    \code{NULL} and selects every parameter. Each quoted string is
    matched to one or more parameter names with the \code{grep}
    function. For example, if the user specifies \code{Parms=c("eta",
      "tau")}, and if the parameter names are beta[1], beta[2], eta[1],
    eta[2], and tau, then all parameters will be selected, because the
    string \code{eta} is within \code{beta}. Since \code{grep} is used,
    string matching uses regular expressions, so beware of
    meta-characters, though these are acceptable: ".", "[", and "]".}
}
\details{\code{PosteriorChecks} is a supplemental function that returns
  a list with two components.}
\value{
  \code{PosteriorChecks} returns an object of class
  \code{posteriorchecks} that is a list with the following components:
  \item{Posterior.Correlation}{
    This is a correlation matrix of the parameters selected with the
    \code{Parms} argument. This component is returned as \code{NA} for
    objects of class \code{"laplace"}.}
  \item{Posterior.Summary}{This is a matrix in which each row is a
    parameter and there are six columns: p(theta > 0), N.Modes,
    Kurtosis, Skewness, Burn-In, and IAT. The first column,
    p(theta > 0), indicates parameter importance by reporting how much
    of the distribution is greater than zero. An important parameter
    distribution will have a result at least as extreme as 0.025 or
    0.975, and an unimportant parameter distribution is centered at
    0.5. This is not the importance of the associated variable relative
    to how well the model fits the data. For variable importance, see
    the \code{\link{Importance}} function. The second column, N.Modes,
    is the number of modes, estimated with the \code{\link{Modes}}
    function. Kurtosis and skewness are useful posterior checks that may
    suggest that a posterior distribution is non-normal or does not fit
    well with a distributional assumption, assuming a distributional
    assumption exists, which it may not. The burn-in is estimated for
    each chain (only for objects of class \code{demonoid}. Lastly, the
    integrated autocorrelation time is estimated with \code{\link{IAT}}.)
    with the \code{\link{burnin}} function.}
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{burnin}},
  \code{\link{IAT}},
  \code{\link{Importance}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{Modes}}.
}
\examples{### See the LaplacesDemon function for an example.}
\keyword{posterior, check}