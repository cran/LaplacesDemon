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
  \item{obj}{This is a vector or matrix of posterior samples.}
  \item{prob}{This is a numeric scalar in the interval (0,1) giving the
    target probability interval.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  For each parameter (or deviance or monitored variable), the HPD is
  estimated from the empirical CDF of the sample as the shortest
  interval for which the difference in the ECDF values of the endpoints
  is the user-specified probability width. This assumes the
  distribution is not severely multi-modal. The HPD contains the
  region of the highest posterior density.

  The HPD is essentially a Bayesian confidence interval. Bayesians
  usually refer to the HPD as a credible interval or probability interval.

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
#HPD(Fit$Posterior1)
}
\seealso{\code{\link{LaplacesDemon}}}