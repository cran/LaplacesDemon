\name{is.demonoid}
\alias{is.demonoid}
\alias{is.laplace}
\title{Demonoid and Laplace}
\description{
  These functions provide a logical test of the class of an object.
}
\usage{
is.demonoid(x)
is.laplace(x)
}
\arguments{
  \item{x}{This is an object that will be tested either for class
    \code{demonoid} or class \code{laplace}.}
}
\details{
  After updating a model with the \code{\link{LaplacesDemon}}
  function, an output object is created. The output object is of class
  \code{demonoid}. Likewise, after passing a model to the
  \code{\link{LaplaceApproximation}} function, an output object is
  created, and it is of class \code{laplace}. The class of these objects
  may be tested either with the \code{is.demonoid} or \code{is.laplace}
  functions, which provide logical tests.
}
\value{
  The \code{is.demonoid} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid}.

  The \code{is.laplace} function returns a logical value indicating
  whether or not the supplied object is of class \code{laplace}.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\keyword{utility}