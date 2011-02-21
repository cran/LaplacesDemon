\name{plot.laplace.ppc}
\alias{plot.laplace.ppc}
\title{Plots of Posterior Predictive Checks}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{laplace.ppc} Density plots are provided.
}
\usage{\method{plot}{laplace.ppc}(x, Rows=NULL, PDF=FALSE, \dots)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{laplace.ppc}.
     } 
     \item{Rows}{
     This optional argument is for a vector of row numbers that 
     specify the distributions associated with the row numbers to plot.
     }
     \item{PDF}{
     This logical argument indicates whether or not the user wants
     Laplace's Demon to save the plots as a .pdf file.
     }
     \item{\dots}{Additional arguments are unused.}
}
\details{
  The density plots show the kernel density of the posterior predictive
  distribution for each selected row of y (all are selected by
  default). A vertical red line indicates the position of the
  observed y along the x-axis. When the vertical red line is close
  to the middle of a normal posterior predictive distribution, then
  there is little discrepancy between y and the posterior predictive
  distribution. When the vertical red line is in the tail of the
  distribution, or outside of the kernel density altogether, then there
  is a large discrepancy between y and the posterior predictive
  distribution. Large discrepancies may be considered outliers, and
  moreover suggest that an improvement in model fit should be
  considered.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{predict.laplace}}}
\examples{### See the LaplaceApproximation function for an example.}
\keyword{plot}