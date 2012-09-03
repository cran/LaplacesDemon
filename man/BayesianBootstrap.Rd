\name{BayesianBootstrap}
\alias{BayesianBootstrap}
\title{The Bayesian Bootstrap}
\description{
  This function samples data for bootstrapping, as per the Bayesian
  Bootstrap of Rubin (1981).
}
\usage{
BayesianBootstrap(X, n=1000, Status=NULL)
}
\arguments{
  \item{X}{This is a vector or matrix of data. When a matrix is
    supplied, sampling is based on the first column.}
  \item{n}{This is the number of bootstrapped samples to draw, and is
    usually equal to the number of rows in \code{X}.}
  \item{Status}{This determines the periodicity of status messages. When
    \code{Status=100}, for example, a status message is displayed every
    100 samples. Otherwise, \code{Status} defaults to \code{NULL}, and
    status messages are not displayed.}
}
\details{
  The term, `bootstrap', comes from the German novel \emph{Adventures of
  Baron Munchausen} by Rudolph Raspe, in which the hero saves himself
  from drowning by pulling on his own bootstraps. The idea of the
  statistical bootstrap is to evaluate properties of an estimator
  through the empirical, rather than theoretical, CDF.
  
  Rubin (1981) introduced the Bayesian bootstrap. In contrast to the
  frequentist bootstrap which simulates the sampling distribution of a
  statistic estimating a parameter, the Bayesian bootstrap simulates the
  posterior distribution.

  The data, \eqn{\textbf{X}}{X}, are assumed to be independent and
  identically distributed (IID), and to be a representative sample of
  the larger (bootstrapped) population. The data are randomly sampled
  with replacement \eqn{n} times. Posterior distributions of interest
  are then estimated given the sampled (bootstrapped) data.

  Just as with the frequentist bootstrap, inappropriate use of the
  Bayesian bootstrap can lead to inappropriate inferences. The Bayesian
  bootstrap violates the Likelihood Principle, because the evaluation of
  a statistic of interest depends on data sets other than the observed
  data set.

  The \code{BayesianBootstrap} function has many uses, including
  imputation (with this variation: \code{\link{ABB}}), validation,
  creating test statistics on the population data given the observed
  data, and more.
}
\value{
  This function returns a \eqn{n \times J}{n x J} matrix of samples,
  where \eqn{n} is specified by the user and \eqn{J} is the number of
  columns in \eqn{\textbf{X}}{X}.
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}

  Rubin, D.B. (1981). "The Bayesian Bootstrap". The Annals of
  Statistics, 9(1), p. 130--134.
}
\seealso{
  \code{\link{ABB}},
  \code{\link{LaplaceApproximation}}, and
  \code{\link{LaplacesDemon}}.
}
\examples{
X <- matrix(runif(100), 10, 10)
X.B <- BayesianBootstrap(X, n=10)
}
\keyword{Utility}
