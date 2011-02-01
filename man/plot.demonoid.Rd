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
\details{
  The plots are arranged in a 3 x 3 matrix.  Each row represents a
  parameter, the deviance, or a monitored variable.  The left column
  displays trace plots, the middle column displays kernel density
  plots, and the right column displays autocorrelation (ACF) plots.

  Trace plots show the thinned history of the chain or Markov chain,
  with its value in the y-axis moving by iteration across the x-axis.
  A chain or Markov chain with good properties does not suggest a trend
  upward or downward as it progresses across the x-axis (it should
  appear stationary), and it should mix well, meaning it should appear
  as though random samples are being taken each time from the same
  target distribution.  Visual inspection of a trace plot cannot
  verify convergence, but apparent non-stationarity or poor mixing can
  certainly suggest non-convergence.

  Kernel density plots depict the marginal posterior distribution.
  There is no distributional assumption about this density.

  Autocorrelation plots show the autocorrelation or serial correlation
  between sampled values at nearby iterations.  Samples with
  autocorrelation do not violate any assumption, but are inefficient
  because they reduce the effective sample size (ESS), and indicate that
  the chain is not mixing well, since each value is influenced by values
  that are previous and nearby.  The x-axis indicates lags with respect
  to samples by iteration, and the y-axis represents autocorrelation.
  The ideal autocorrelation plot shows perfect correlation at zero lag,
  and quickly falls to zero autocorrelation for all other lags.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
