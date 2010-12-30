\name{predict.demonoid}
\alias{predict.demonoid}
\title{Posterior Predictive Checks}
\description{
  This may be used to perform posterior predictive checks, wherein
  replicates of y or y[rep] are predicted, given an object of class
  \code{demonoid}, the model specification, and data.
}
\usage{\method{predict}{demonoid}( object, Log.Posterior, Data, \dots)}
\arguments{
     \item{object}{An object of class \code{demonoid} is required.}
     \item{Log.Posterior}{The model specification function is required.}
     \item{Data}{A data set in a list is required.}
     \item{\dots}{Additional arguments are unused.}
}
\details{
  This function performs posterior predictive checks.  The purpose of a
  posterior predictive check is to assess how well (or poorly) the model
  fits the data, or to assess discrepancies between the model and the
  data.  Several methods exist for posterior predictive checks, but 
  currently this function only predicts y[rep] and reports concordance 
  between y and y[rep].

  At the moment, y must be univariate, and this function does not
  predict y[new], or unobserved instances of y.
}
\value{
     This function returns a list with the following items:
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
          This is a n x 6 matrix, where n is the number of records of
          y and there are 6 columns, as follows: Mean, SD, LB (the 2.5\%
	  quantile), Median, UB (the 97.5\% quantile), and PPPV, the
	  posterior predictive p-value, which is Pr(y[rep] >= y), as in
	  Gelman, Meng, and Stern (1996).
          }
     \item{yrep}{
          This is a n x s matrix, where n is the number of records of
          y and s is the number of posterior samples.
          }
}
\references{
     Gelman, A., Meng, X.L., and Stern, H. (1996). Posterior Predictive
     Assessment of Model Fitness via Realized Discrepancies. Statistica
     Sinica, 6, p. 733-807.

     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
