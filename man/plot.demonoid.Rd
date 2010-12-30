\name{plot.demonoid}
\alias{plot.demonoid}
\title{Plot samples from the output of \code{LaplacesDemon}}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{demonoid}  Plots include a trace plot, density plot, and
  autocorrelation or ACF plot.
}
\usage{\method{plot}{demonoid}( x, BurnIn=1, Data, PDF=FALSE, Parms, \dots)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{demonoid}.
     } 
     \item{BurnIn}{
     This argument requires a positive integer that indicates a 
     number of iterations to discard as burn-in for the purposes of
     plotting.
     }
     \item{Data}{
     This required argument must receive the list of data that
     was supplied to \code{LaplacesDemon} to create the object of
     class \code{demonoid}.
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
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
