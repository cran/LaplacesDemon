\name{summary.demonoid.ppc}
\alias{summary.demonoid.ppc}
\title{Posterior Predictive Check Summary}
\description{
  This may be used to summarize either new, unobserved instances of y
  (called y[new]) or replicates of y (called y[rep]). Either y[new] or
  y[rep] is summarized, depending on \code{predict.demonoid}.
}
\usage{\method{summary}{demonoid.ppc}(object, Rows, Discrep, d, \dots)}
\arguments{
     \item{object}{An object of class \code{demonoid.ppc} is required.}
     \item{Rows}{An optional vector of row numbers, for example
          \code{c(1:10)}. All rows will be estimated, but only these
          rows will appear in the summary.}
     \item{Discrep}{A character string indicating a discrepancy test.
          \code{Discrep} defaults to \code{NULL}. Valid character
          strings are: "max(yhat[i,]) > max(y)",
	  "min(yhat[i,]) < min(y)", "round(yhat[i,]) = d".}
     \item{d}{This is an optional integer to be used with the
	  \code{Discrep} argument above, and it defaults to \code{d=0}.}
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

  A discrepancy measure is an approach to studying discrepancies between
  the model and data (Gelman et al., 1996).
  
  In the future, other discrepancy measures will be included.
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
     \item{Discrepancy.Statistic}{
          This is only reported if the \code{Discrep} argument receives a
          valid discrepancy measure as listed above. The \code{Discrep}
	  applies to each record of y, and the
	  \code{Discrepancy.Statistic} reports the results of the
	  discrepancy measure on the entire data set. For example, if
	  \code{Discrep="min(yhat[i,] < min(y)"}, then the overall
	  result is the proportion of records in which the minimum
	  sample of yhat was less than the overall minimum y. This is
	  Pr(min(yhat[i,]) < min(y) | y, Theta), where Theta is the
	  parameter set.}
     \item{Summary}{
          This is a n x 8 matrix, where n is the number of records of
          y and there are 8 columns, as follows: y, Mean, SD, LB (the
	  2.5\% quantile), Median, UB (the 97.5\% quantile), PQ (the
	  predictive quantile, which is Pr(y[rep] >= y)), and Test,
	  which shows the record-level result of a test, if specified.
          }
}
\references{
Gelfand, A. (1996). "Model Determination Using Sampling Based Methods". 
In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in 
Markov Chain Monte Carlo in Practice. Chapman \& Hall: Boca 
Raton, FL.

Gelman, A., Meng, X.L., and Stern H. (1996). "Posterior Predictive
Assessment of Model Fitness via Realized Discrepancies". Statistica
Sinica, 6, p. 733--807.

Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{predict.demonoid}}}
\examples{### See the LaplacesDemon function for an example.}
\keyword{summary}