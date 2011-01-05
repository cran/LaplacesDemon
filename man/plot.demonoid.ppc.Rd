\name{plot.demonoid.ppc}
\alias{plot.demonoid.ppc}
\title{Plots of Posterior Predictive Checks}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{demonoid.ppc}  Plots include a density plot.
}
\usage{\method{plot}{demonoid.ppc}( x, Rows=NULL, PDF=FALSE, \dots)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{demonoid.ppc}.
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
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{predict.demonoid}}}
