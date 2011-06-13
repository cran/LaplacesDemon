\name{caterpillar.plot}
\alias{caterpillar.plot}
\title{Caterpillar Plot}
\description{
  This may be used to produce a caterpillar plot of posterior samples
  from an object either of class \code{demonoid} or class
  \code{laplace}. A caterpillar plot is a horizontal plot of 3
  quantiles of selected distributions.
}
\usage{caterpillar.plot(x, Parms=NULL, Title=NULL)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{demonoid} or
     class \code{laplace}. For an object of class \code{demonoid}, the
     distributions of the stationary posterior summary (\code{Summary2})
     will be attempted first, and if missing, then the parameters of all
     posterior samples (\code{Summary1}) will be plotted. For an object
     of class \code{laplace}, the distributions in the posterior
     summary, \code{Summary}, are plotted.
     }
     \item{Parms}{
       This argument accepts a concatenation of rows of the posterior
       summary tables to be plotted. For example, \code{Parms=c(1,3,5)}
       will plot the first, third, and fifth parameters listed in the
       posterior summary table. This argument defaults to \code{NULL},
       which will plot all parameters, including the deviance and
       monitored variables.
     }
     \item{Title}{
       This argument accepts a title for the plot.
     }
   }
\details{
  Caterpillar plots are popular plots in Bayesian inference for
  summarizing the quantiles of posterior samples. A caterpillar plot is
  similar to a horizontal boxplot, though without quartiles, making it
  easier for the user to study more distributions in a single plot. The
  following quantiles are plotted as a line for each parameter: 0.025 and
  0.975. The median appears as a black dot. A vertical, gray line
  is included at zero.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{#An example is provided in the LaplacesDemon function.}
