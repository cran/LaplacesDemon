\name{summary.demonoid.ppc}
\alias{summary.demonoid.ppc}
\title{Posterior Predictive Check Summary}
\description{
  This may be used to summarize either new, unobserved instances of y
  (called y[new]) or replicates of y (called y[rep]). Either y[new] or
  y[rep] is summarized, depending on \code{predict.demonoid}.
}
\usage{\method{summary}{demonoid.ppc}( object, Rows, \dots)}
\arguments{
     \item{object}{An object of class \code{demonoid.ppc} is required.}
     \item{Rows}{An optional vector of row numbers, for example
          \code{c(1:10)}. All rows will be estimated, but only these
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
  the model and the data. When y is known, this function estimates the
  predictive concordance between y and y[rep] as per Gelfand (1996), and
  the predictive quantile (PQ), which is the individual-level metric for
  outlier detection used to calculate Gelfand's predictive concordance.

  At the moment, y must be univariate and continuous.

  Numerous other posterior predictive checks exist. A future goal is to
  include the Conditional Predictive Ordinate (CPO).
}
\value{
     This function returns a list with the following components:
     \item{Concordance}{
          This is the percentage of the records of y that are within the 
          95\% probability interval of y[rep]. A probability interval
	  is also called a credible interval. Gelfand's suggested goal
	  is to achieve 95\% predictive concordance. Lower percentages
	  indicate too many outliers and a poor fit of the model to the
	  data, and higher percentages may suggest overfitting.
          }
     \item{Summary}{
          This is a n x 7 matrix, where n is the number of records of
          y and there are 7 columns, as follows: y, Mean, SD, LB (the
	  2.5\% quantile), Median, UB (the 97.5\% quantile), and
	  PQ, the predictive quantile, which is Pr(y[rep] >= y).
          }
}
\references{
Gelfand, A. (1996). "Model Determination Using Sampling Based Methods". 
In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in 
Markov Chain Monte Carlo in Practice. Chapman \& Hall: Boca 
Raton, FL.

Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{predict.demonoid}}}
