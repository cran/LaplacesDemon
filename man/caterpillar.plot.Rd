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
     This argument accepts a vector of quoted strings to be matched for
     selecting parameters for plotting. This argument defaults to
     \code{NULL} and selects every parameter for plotting. Each quoted
     string is matched to one or more parameter names with the
     \code{grep} function. For example, if the user specifies
     \code{Parms=c("eta", "tau")}, and if the parameter names
     are beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
     will be selected, because the string \code{eta} is within
     \code{beta}. Since \code{grep} is used, string matching uses
     regular expressions, so beware of meta-characters, though these are
     acceptable: ".", "[", and "]".
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
