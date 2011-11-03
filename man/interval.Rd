\name{interval}
\alias{interval}
\title{Constrain to Interval}
\description{
  This function constrains the value(s) of a scalar, vector, or matrix
  to a specified interval, \eqn{[a,b]}. In Bayesian inference, it is
  often useful to truncate a prior distribution to an interval, such as
  \eqn{p(\theta) \in [a,b]}{p(theta) in [a,b]}. The \code{interval}
  function is often used in conjunction with the \code{\link{dtrunc}}
  function to truncate a probability distribution. While
  \code{\link{dtrunc}} prevents assigning density outside of its
  interval and re-estimates density within the interval, the
  \code{interval} function is used to prevent the parameter from moving
  outside of the interval in the first place.

  After the parameter is constrained to an interval in
  \code{\link{LaplaceApproximation}} or \code{\link{LaplacesDemon}}, the
  constrained parameter should be updated back into the \code{parm}
  vector, so the algorithm knows it has been constrained. In
  \code{\link{LaplaceApproximation}}, however, note that
  \code{Method="CG"} will be unaware from iteration to iteration that
  the user constrained a parameter to an interval, so the much slower
  algorithm, \code{Method="AGA"}, must be used instead when constraining
  intervals.
}
\usage{
interval(x, a=-Inf, b=Inf)
}
\arguments{
  \item{x}{This is a scalar, vector, or matrix, and its elements will be
    constrained to the interval [\code{a},\code{b}].}
  \item{a}{This optional argument allows the specification of the lower
    bound of the interval, and defaults to \code{-Inf}.}
  \item{b}{This optional argument allows the specification of the upper
    bound of the interval, and defaults to \code{Inf}.}
}
\details{
  For \code{\link{LaplacesDemon}}, the MCMC algorithm makes a proposal
  that may or may not be accepted. If it is accepted, then the accepted
  values may be constrained in the model specification, and returned
  to the algorithm with \code{parm}. This is an inconsistency between
  what the algorithm accepted, and what was used. Perhaps the
  constrained values would have been rejected. After running hundreds of
  different types of models in Laplace's Demon with constrained
  parameters, it seems this inconsistency does not pose any problems,
  but it is important to note it.
}
\value{
  The \code{interval} function returns a scalar, vector, or matrix in
  accord with its argument, \code{x}. Each element is constrained to
  the interval [\code{a},\code{b}].
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{dtrunc}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\examples{
#See the Examples vignette for numerous examples.
library(LaplacesDemon)
x <- 2
interval(x,0,1)
X <- matrix(runif(25,-2,2),5,5)
interval(X,-1,1)
}
\keyword{utility}