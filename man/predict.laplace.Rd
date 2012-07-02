\name{predict.laplace}
\alias{predict.laplace}
\title{Posterior Predictive Checks}
\description{
  This may be used to predict either new, unobserved instances of
  \eqn{\textbf{y}}{y} (called \eqn{\textbf{y}^{new}}{y[new]}) or
  replicates of \eqn{\textbf{y}}{y} (called
  \eqn{\textbf{y}^{rep}}{y[rep]}), and then perform posterior
  predictive checks. Either \eqn{\textbf{y}^{new}}{y[new]} or
  \eqn{\textbf{y}^{rep}}{y[rep]} is predicted given an object of
  class \code{laplace}, the model specification, and data. This function
  requires that posterior samples were produced with
  \code{\link{LaplaceApproximation}}.
}
\usage{\method{predict}{laplace}(object, Model, Data, \dots)}
\arguments{
  \item{object}{An object of class \code{laplace} is required.}
  \item{Model}{The model specification function is required.}
  \item{Data}{A data set in a list is required. The dependent
     variable is required to be named either \code{y} or \code{Y}.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  Since Laplace Approximation characterizes marginal posterior
  distributions with modes and variances, and posterior predictive
  checks involve samples, the \code{predict.laplace} function requires
  the use of independent samples of the marginal posterior
  distributions, provided by \code{\link{LaplaceApproximation}} when
  \code{sir=TRUE}.

  The samples of the marginal posterior distributions of the target
  distributions (the parameters) are passed along with the data to the
  \code{Model} specification and used to draw samples from the deviance
  and monitored variables. At the same time, the fourth component in the
  returned list, which is labeled \code{yhat}, is a vector of
  expectations of \eqn{\textbf{y}}{y}, given the samples, model
  specification, and data. To predict \eqn{\textbf{y}^{rep}}{y[rep]},
  simply supply the data set used to estimate the model. To predict
  \eqn{\textbf{y}^{new}}{y[new]}, supply a new data set instead (though
  for some model specifications, this cannot be done, and
  \eqn{\textbf{y}_{new}}{y[new]} must be specified in the \code{Model}
  function). If the new data set does not have \eqn{\textbf{y}}{y}, then
  create \code{y} in the list and set it equal to something sensible,
  such as \code{mean(y)} from the original data set.

  The variable \code{y} must be a vector. If instead it is matrix
  \code{Y}, then it will be converted to vector \code{y}. The vectorized
  length of \code{y} or \code{Y} must be equal to the vectorized length
  of \code{yhat}, the fourth component of the returned list of the
  \code{Model} function.
}
\value{
  This function returns an object of class \code{laplace.ppc} (where
  ``ppc'' stands for posterior predictive checks). The returned object
  is a list with the following components:
  \item{y}{
    This stores \eqn{\textbf{y}}{y}, the dependent variable.}
  \item{yhat}{
    This is a \eqn{N \times S}{N x S} matrix, where \eqn{N} is the
    number of records of \eqn{\textbf{y}}{y} and \eqn{S} is the number of 
    posterior samples.}
  \item{deviance}{
    This is a vector of length \eqn{S}, where \eqn{S} is the number of
    independent posterior samples. Samples are obtained with the
    sampling importance resampling algorithm, \code{\link{SIR}}.}
  \item{monitor}{
    This is a \eqn{N \times S}{N x S} matrix, where \eqn{N} is the
    number of monitored variables and \eqn{S} is the number of independent
    posterior samples. Samples are obtained with the sampling importance
    resampling algorithm, \code{\link{SIR}}.}
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{SIR}}.
}
\keyword{Posterior Predictive Checks}
\keyword{Predict}
