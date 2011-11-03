\name{Consort}
\alias{Consort}
\title{Consort with Laplace's Demon}
\description{
  This may be used to consort with Laplace's Demon regarding an object
  of class \code{demonoid}. Laplace's Demon will offer suggestions.
}
\usage{Consort(object)}
\arguments{
     \item{object}{
     This required argument is an object of class \code{demonoid}.
     }
}
\details{
  First, \code{Consort} calls \code{print.demonoid}, which prints most of
  the components to the screen from the supplied object of class
  \code{demonoid}.

  Second, Laplace's Demon considers a combination of five conditions
  when making the largest part of its suggestion. These conditions are:
  the algorithm, acceptance rate, MCSE, ESS, and stationarity. There
  are 48 combinations of these five conditions, though several
  combinations yield the same suggestion. Other things are considered
  as well, such as the recommended thinning value is used to suggest a
  new number of iterations. The \code{Periodicity} argument is
  suggested to be set equal to the value of \code{Rec.Thinning}.

  Appeasement Conditions
  \itemize{
  \item Algorithm: The final algorithm must be non-adaptive, so that the
    Markov property holds.
  \item Acceptance Rate: The acceptance rate is considered satisfactory
    if it is within the interval [15\%,50\%].
  \item MCSE: The Monte Carlo Standard Error (MCSE) is considered
    satisfactory for each target distribution if it is less than 6.27\%
    of the standard deviation of the target distribution. This allows
    the true mean to be within 5\% of the area under a Gaussian
    distribution around the estimated mean.
  \item ESS: The effective sample size (ESS) is considered satisfactory
    for each target distribution if it is at least 100, which is usually
    enough to describe 95\% probability intervals. The \code{\link{ESS}}
    function is used.
  \item Stationarity: Each target distribution is considered
    satisfactory if it is estimated to be stationary with the
    \code{\link{Geweke.Diagnostic}} function.
}

  If all five conditions are satisfactory, then Laplace's Demon is
  appeased. Otherwise, Laplace's Demon will suggest and supply R
  code that is ready to be copy/pasted and executed.

  The Demonic Suggestion is usually very helpful, but should not be
  followed blindly. Do not let it replace critical thinking.
}
\references{
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{ESS}},
  \code{\link{Geweke.Diagnostic}},
  \code{\link{LaplacesDemon}}.
}
