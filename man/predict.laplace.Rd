\name{predict.laplace}
\alias{predict.laplace}
\title{Posterior Predictive Checks}
\description{
  This may be used to predict either new, unobserved instances of y
  (called y[new]) or replicates of y (called y[rep]), and then perform
  posterior predictive checks. Either y[new] or y[rep] is predicted
  given an object of class \code{laplace}, the model specification, and
  data. This function also samples the deviance and monitored variables.
}
\usage{\method{predict}{laplace}(object, Model, Data, Samples, \dots)}
\arguments{
     \item{object}{An object of class \code{laplace} is required.}
     \item{Model}{The model specification function is required.}
     \item{Data}{A data set in a list is required. The dependent
          variable is required to be named either \code{y} or \code{Y}.}
     \item{Samples}{This is the number of posterior samples to be created.}
     \item{\dots}{Additional arguments are unused.}
}
\details{
  This function passes each iteration of marginal posterior samples
  along with data to \code{Model}, where the fourth component in the
  return list is labeled \code{yhat}, and is a vector of expectations of
  y, given the samples, model specification, and data. To predict
  y[rep], simply supply the data set used to estimate the model. To
  predict y[new], supply a new data set instead. If the new data set
  does not have y, then create \code{y} in the list and set it equal to
  something sensible, such as \code{mean(y)} from the original data
  set.

  The variable \code{y} must be a vector. If instead it is matrix
  \code{Y}, then it will be converted to vector \code{y}. The vectorized
  length of \code{y} or \code{Y} must be equal to the vectorized length
  of \code{yhat}, the fourth component of the return list of the
  \code{Model} function.

  This function also samples the deviance and monitored variables.
}
\value{
     This function returns an object of class \code{laplace.ppc} (where
     ppc stands for posterior predictive checks). The returned object is
     a list with the following components:
     \item{y}{
          This stores y, the dependent variable.
          }
     \item{yhat}{
          This is a n x s matrix, where n is the number of records of
          y and s is the number of posterior samples.
          }
     \item{deviance}{
          This is a vector of length s, where s is the number of
	  posterior samples.
          }
     \item{monitor}{
          This is a n x s matrix, where n is the number of monitored
	  variables and s is the number of posterior samples.
          }
}
\references{
     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{LaplaceApproximation}}}
\keyword{predict}
