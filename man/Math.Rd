\name{Math}
\alias{logadd}
\alias{partial}
\title{Math Utility Functions}
\description{
  These are utility functions for math.
}
\usage{
logadd(x, add=TRUE)
partial(Model, parm, Data, Interval=1e-6)
}
\arguments{
  \item{x}{This is a numeric vector.}
  \item{add}{Logical. This defaults to \code{TRUE}, in which case
    \eqn{\log(x+y)}{log(x+y)} is performed. Otherwise,
    \eqn{\log(x-y)}{log(x-y)} is performed.}
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{parm}{This is a vector parameters.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplacesDemon}}.}
  \item{Interval}{This is the interval of numeric differencing.}
}
\details{
  The \code{logadd} function performs addition (or subtraction) when the
  terms are logarithmic. The equations are:

  \deqn{\log(x+y) = \log(x) + \log(1 + \exp(\log(y) - \log(x)))}{log(x+y) =
    log(x) + log(1 + exp(log(y) - log(x)))}
  \deqn{\log(x-y) = \log(x) + \log(1 - \exp(\log(y) - \log(x)))}{log(x+y) =
    log(x) + log(1 - exp(log(y) - log(x)))}
  
  The \code{partial} function estimates partial derivatives of
  parameters in a model specification with data, using
  finite-differencing. In calculus, a partial derivative of a function
  of several variables is its derivative with respect to one of those
  variables, with the others held constant. In this function, first
  partial derivatives are approximated with numerical differencing, by
  measuring the change in the function with respect to a small change in
  each parameter. A related function is \code{Hessian}, which returns a
  matrix of second-order partial derivatives of the model specification
  function with respect to its parameters. This function is not intended
  to be called by the user, but is used by other functions. Automatic
  differentiation (AD) is an improvement over finite-differencing, but
  is currently unavailable here.
}
\value{
  \code{logadd} returns the result of \eqn{\log(x+y)}{log(x+y)} or
  \eqn{\log(x-y)}{log(x-y)}.
  
  \code{partial} returns a vector of partial derivatives.
}
\author{Statisticat, LLC. \email{statisticat@gmail.com}}
\seealso{
  \code{\link{Hessian}},
  \code{\link{LaplaceApproximation}}, and
  \code{\link{LaplacesDemon}}.
}
\keyword{Math}