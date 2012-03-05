\name{partial}
\alias{partial}
\title{First Partial Derivative}
\description{
This function returns a vector of first partial derivatives of a model
specification function with respect to its parameters.
}
\usage{
partial(Model, parm, Data, Interval=1e-6)
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{parm}{This is a vector parameters.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplacesDemon}}.}
  \item{Interval}{This is the interval of numeric differencing.}
}
\details{
  In calculus, a partial derivative of a function of several variables
  is its derivative with respect to one of those variables, with the
  others held constant. In this function, first partial derivatives are
  approximated with numerical differencing, by measuring the change in
  the function with respect to a small change in each parameter. A
  related function is \code{Hessian}, which returns a matrix of
  second-order partial derivatives of the model specification function
  with respect to its parameters.

  This function is not intended to be called by the user, but is used by
  other functions.
}
\value{
  \code{partial} returns a vector of partial derivatives.
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{Hessian}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\keyword{math}