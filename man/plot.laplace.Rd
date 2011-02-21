\name{plot.laplace}
\alias{plot.laplace}
\title{Plot the output of \code{LaplaceApproximation}}
\description{
  This may be used to plot, or save plots of, the iterated history of
  the parameters in an object of class \code{laplace}.
}
\usage{\method{plot}{laplace}(x, Data, PDF=FALSE, Parms, \dots)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{laplace}.
     }
     \item{Data}{
     This required argument must receive the list of data that
     was supplied to \code{LaplaceApproximation} to create the object
     of class \code{laplace}.
     }
     \item{PDF}{
     This logical argument indicates whether or not the user wants
     Laplace's Demon to save the plots as a .pdf file.
     }
     \item{Parms}{
     This argument accepts a positive integer that indicates how
     many parameters should be plotted, in case there are numerous
     parameters and only a quick glance is desired.
     }
     \item{\dots}{Additional arguments are unused.}
   }
\details{
  The plots are arranged in a 3 x 3 matrix. The monitored variables are
  not plotted. The purpose of these plots is to show how the value of
  each parameter and the deviance changed by iteration as the
  \code{LaplaceApproximation} attempted to maximize the logarithm of the
  unnormalized joint posterior density.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplaceApproximation}}}
\examples{### See the LaplaceApproximation function for an example.}
\keyword{plot}