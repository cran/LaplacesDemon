\name{Consort}
\alias{Consort}
\title{Consort with Laplace's Demon}
\description{
  This may be used to consort with Laplace's Demon regarding an object
  of class \code{demonoid}. Laplace's Demon will offer suggestions.
}
\usage{Consort( object)}
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

  If all five conditions are satisfactory, then Laplace's Demon is
  appeased. Otherwise, Laplace's Demon will suggest and supply R
  code that is ready to be copy/pasted and executed.

  The Demonic Suggestion is usually very helpful, but should not be
  followed blindly. Do not let it replace critical thinking.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
