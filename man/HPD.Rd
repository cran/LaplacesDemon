\name{HPD}
\alias{HPD}
\title{Highest Posterior Density}
\description{
  This function estimates the interval of the Highest Posterior Density (HPD)
  of posterior samples.
}
\usage{
HPD(obj, prob=0.95, \dots)
}
\arguments{
  \item{obj}{This can be either a vector or matrix of posterior samples,
    or an object of class \code{demonoid}. If it is an object of class
    \code{demonoid}, then it will use only stationary posterior samples
    and monitored target distributions (automatically discarding the
    burn-in; if stationarity does not exist, then it will use all
    samples).}
  \item{prob}{This is a numeric scalar in the interval (0,1) giving the
    target probability interval, and defaults to 0.95, representing a
    95\% probability interval.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The \code{HPD} is one form of a Bayesian confidence interval. Bayesians
  usually refer to the \code{HPD} as a credible interval or probability
  interval. The \code{HPD} is identical to a quantile-based probability
  interval when the posterior distribution is uni-modal and
  symmetric. As an example, imagine an exponential posterior
  distribution. A quantile-based probability interval would report the
  highest density region near zero to be outside of its interval. In
  contrast, the \code{HPD} is recommended for such skewed posterior
  distributions. However, posterior distributions with multiple modes
  probably should not be characterized with any probability
  interval. The presence of multiple modes may be checked with the
  \code{\link{is.multimodal}} function.
  
  For each parameter (or deviance or monitored variable), the HPD is
  estimated from the empirical CDF of the sample as the shortest
  interval for which the difference in the ECDF values of the end-points
  is the user-specified probability width. This assumes the
  distribution is not severely multi-modal. The HPD contains the
  region of the highest posterior density.
}
\value{
  A matrix is returned with rows corresponding to the parameters (or
  deviance or monitored variables), and columns \code{"Lower"} and
  \code{"Upper"}. The attribute \code{"Probability"} is the
  user-selected probability width.
}
\examples{
#First, update the model with the LaplacesDemon function.
#Then
#HPD(Fit)
}
\seealso{
  \code{\link{is.multimodal}},
  \code{\link{LaplacesDemon}}.}