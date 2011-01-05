\name{summary.demonoid.ppc}
\alias{summary.demonoid.ppc}
\title{Posterior Predictive Check Summary}
\description{
  This may be used to summarize either new, unobserved instances of y
  (called y[new]) or replicates of y (called y[rep]).  Either y[new] or
  y[rep] is summarized, depending on \code{predict.demonoid}.
}
\usage{\method{summary}{demonoid.ppc}( object, Rows, \dots)}
\arguments{
     \item{object}{An object of class \code{demonoid.ppc} is required.}
     \item{Rows}{An optional vector of row numbers, for example
          \code{c(1:10)}.  All rows will be estimated, but only these
          rows will appear in the summary.}
     \item{\dots}{Additional arguments are unused.}
}
\details{
  This function summarizes an object of class \code{demonoid.ppc}, which
  consists of posterior predictive checks on either y[new] or y[rep],
  depending respectively on whether unobserved instances of y or the
  model sample of y was used in the \code{predict.demonoid} function.
  
  The purpose of a posterior predictive check is to assess how well (or
  poorly) the model fits the data, or to assess discrepancies between
  the model and the data.  When y is known, this function estimates a
  chi-square, the concordance between y and y[rep], and the posterior
  predictive p-value (p.value), as in Gelman, Meng, and Stern (1996).

  At the moment, y must be univariate and continuous.

  Other popular posterior predictive checks that would be nice to
  include in the future are the Posterior Predictive Ordinate (PPO) and
  Conditional Predictive Ordinate (CPO).
}
\value{
     This function returns a list with the following items:
     \item{Chi.Square}{
          This is a chi-square to represent model fit (Gelman, Carlin,
          Stern, and Rubin, 2004).
          }
     \item{Concordance}{
          This is the percentage of the records of y that are within the 
          95\% probability interval of y[rep].  A probability interval
	  is also called a credible interval.  A model that fits the
	  data well has few instances of y outside of the 95\% probability
	  interval of y[rep], and each instance of y outside of this
	  interval may be considered an outlier in fit.  If a pattern
	  exists, then it should be studied and the model refined.
          }
     \item{Summary}{
          This is a n x 7 matrix, where n is the number of records of
          y and there are 7 columns, as follows: y, Mean, SD, LB (the
	  2.5\% quantile), Median, UB (the 97.5\% quantile), and
	  p.value, the posterior predictive p-value, which is
	  Pr(y[rep] >= y).
          }
}
\references{
     Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). Bayesian
     Data Analysis, Texts in Statistical Science, 2nd ed., Chapman and
     Hall, London.
     
     Gelman, A., Meng, X.L., and Stern, H. (1996). Posterior Predictive
     Assessment of Model Fitness via Realized Discrepancies. Statistica
     Sinica, 6, p. 733-807.

     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{predict.demonoid}}}
