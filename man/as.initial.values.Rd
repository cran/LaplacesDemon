\name{as.initial.values}
\alias{as.initial.values}
\title{Initial Values}
\description{
This function returns the most recent posterior samples from an object
of class \code{demonoid}, or the posterior modes of an object of class
\code{laplace}, to be the initial values for future updates.
}
\usage{
as.initial.values(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid} or class
    \code{laplace}.}
}
\value{
     \code{as.initial.values} returns the most recent posterior samples
     from an object of class \code{demonoid}, or the posterior modes of
     an object of class \code{laplace}, to be the initial values for
     future updates.
}
\references{
     Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\keyword{initial values}