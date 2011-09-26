\name{summary.laplace.ppc}
\alias{summary.laplace.ppc}
\title{Posterior Predictive Check Summary}
\description{
  This may be used to summarize either new, unobserved instances of y
  (called y[new]) or replicates of y (called y[rep]), as well as the
  deviance and monitored variables. These variables are summarized,
  depending on \code{predict.laplace}.
}
\usage{\method{summary}{laplace.ppc}(object, Categorical, Rows, Discrep,
     d, Quiet, \dots)}
\arguments{
     \item{object}{An object of class \code{laplace.ppc} is required.}
     \item{Categorical}{Logical. If \code{TRUE}, then \code{y} and
          \code{yhat} are considered to be categorical (such as y=0 or
           y=1),  rather than continuous.}
     \item{Rows}{An optional vector of row numbers, for example
          \code{c(1:10)}. All rows will be estimated, but only these
          rows will appear in the summary.}
     \item{Discrep}{A character string indicating a discrepancy test.
          \code{Discrep} defaults to \code{NULL}. Valid character
          strings when \code{y} is continuous are: "Chi-Square",
	  "Kurtosis", "Skewness", "max(yhat[i,]) > max(y)",
	  "min(yhat[i,]) < min(y)", "round(yhat[i,]) = d". Valid
	  character strings when \code{y} is categorical are:
	  "p(yhat[i,] != y[i])". Kurtosis and skewness are not
	  discrepancies, but are included here for convenience.}
     \item{d}{This is an optional integer to be used with the
	  \code{Discrep} argument above, and it defaults to \code{d=0}.}
     \item{Quiet}{This logical argument defaults to \code{FALSE} and
          will print results to the console. When \code{TRUE},
          results are not printed.}
     \item{\dots}{Additional arguments are unused.}
}
\details{
  This function summarizes an object of class \code{laplace.ppc}, which
  consists of posterior predictive checks on either y[new] or y[rep],
  depending respectively on whether unobserved instances of y or the
  model sample of y was used in the \code{predict.laplace} function. The
  deviance and monitored variables are also summarized, though these are
  not posterior predictive checks.
  
  The purpose of a posterior predictive check is to assess how well (or
  poorly) the model fits the data, or to assess discrepancies between
  the model and the data.

  When y is continuous and known, this function estimates the predictive
  concordance between y and y[rep] as per Gelfand (1996), and the
  predictive quantile (PQ), which is the individual-level metric for
  outlier detection used to calculate Gelfand's predictive concordance.

  When y is categorical and known, this function estimates the
  record-level lift, which is p(yhat[i,] = y[i]) / [p(y = j) / n], or
  the number of correctly predicted samples over the rate of that
  category of y in vector y.

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
	  Concordance occurs only when \code{y} is continuous.
          }
     \item{Mean Lift}{
	  This is the mean of the record-level lifts, and occurs only
	  when \code{y} is categorical.}
     \item{DIC (Dbar)}{
          This is the mean deviance, and is a component of the Deviance
          Information Criterion.
          }
     \item{DIC (pD)}{
          This is the effective number of parameters, and is a component
          of the Deviance Information Criterion. \code{pD} is calculated
	  as \code{var(Deviance) / 2}.
          }
     \item{DIC (DIC)}{
          This is the Deviance Information Criterion, a measure of model
          fit. It is the sum of Dbar and pD.
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
	\item{L-criterion}{
	  The L-criterion (Laud and Ibrahim, 1995) was developed for
	  model and variable selection. It is a sum of two components:
	  one involves the predictive variance and the other includes
	  the accuracy of the means of the predictive distribution. The
	  L-criterion measures model performance with a combination of
	  how close its predictions are to the observed data and
	  variability of the predictions. Better models have smaller
	  values of \code{L}. \code{L} is measured in the same units as
	  the response variable, and measures how close the data vector
	  \code{y} is to the predictive distribution. In addition to the
	  value of \code{L}, there is a value for \code{S.L}, which is
	  the calibration number of \code{L}, and is useful in
	  determining how much of a decrease is necessary between models
	  to be noteworthy.}
     \item{Monitor}{
          This is a n x 5 matrix, where n is the number of monitored 
          variables and there are 5 columns, as follows: Mean, SD, LB (the
	  2.5\% quantile), Median, and UB (the 97.5\% quantile).
          }
     \item{Summary}{
	  When \code{y} is continuous, this is a n x 8 matrix, where n
          is the number of records of \code{y} and there are 8 columns,
	  as follows: y, Mean, SD, LB (the 2.5\% quantile), Median, UB
	  (the 97.5\% quantile), PQ (the predictive quantile, which is
	  Pr(y[rep] >= y)), and Test, which shows the record-level
	  result of a test, if specified. When \code{y} is categorical,
	  this matrix has a number of columns equal to the number of
	  categories of \code{y} plus 3, also including \code{y},
	  \code{Lift}, and \code{Discrep}.
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

Laud, P.W. and Ibrahim, J.G. (1995). "Posterior Model Selection".
Journal of the Royal Statistical Society, B 57, p. 247--262.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{predict.laplace}}}
\examples{### See the LaplaceApproximation function for an example.}
\keyword{summary}