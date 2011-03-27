\name{interval}
\alias{interval}
\title{Constrain to Interval}
\description{
This function constrains the value(s) of a scalar, vector, or matrix to
an interval.
}
\usage{
interval(x, a=-Inf, b=Inf)
}
\arguments{
    \item{x}{This is a scalar, vector, or matrix, and its elements will
      be constrained to the interval [\code{a},\code{b}].}
    \item{a}{This optional argument allows the specification of the
      lower bound of the interval, and defaults to \code{-Inf}.}
    \item{b}{This optional argument allows the specification of the
      upper bound of the interval, and defaults to \code{Inf}.}
}
\value{
     \code{interval} returns a scalar, vector, or matrix in accord with
     its argument \code{x}. Each element is constrained to the interval [\code{a},\code{b}].
}
\examples{

library(LaplacesDemon)
x <- 2
interval(x,0,1)
X <- matrix(runif(25,-2,2),5,5)
interval(X,-1,1)
}
\keyword{utility}