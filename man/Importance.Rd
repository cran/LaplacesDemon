\name{Importance}
\alias{Importance}
\title{Variable Importance}
\description{
The \code{Importance} function considers variable importance (or
predictor importance) to be the effect that the variable has on
replicates \eqn{\textbf{y}^{rep}}{y^rep} (or
\eqn{\textbf{Y}^{rep}}{Y^rep}) when the variable is removed
from the model by setting it equal to zero. Here, variable
importance is considered in terms of the comparison of posterior
predictive checks. This may be considered to be a form of sensitivity
analysis, and can be useful for model revision, variable selection, and
model interpretation.

Currently, this function only tests the variable importance of design
matrix \eqn{\textbf{X}}{X}.
}
\usage{
Importance(object, Model, Data, Categorical=FALSE, Discrep, d=0)
}
\arguments{
  \item{object}{An object of class \code{demonoid} or \code{laplace} is
    required.}
  \item{Model}{The model specification function is required.}
  \item{Data}{A data set in a list is required. The dependent variable
    is required to be named either \code{y} or \code{Y}. The
    \code{Importance} function will sequentially remove each column
    vector in \code{X}, so \code{X} is required to be in data set
    \code{Data}.}
  \item{Categorical}{Logical. If \code{TRUE}, then \code{y} and
    \code{yhat} are considered to be categorical (such as y=0 or y=1),
    rather than continuous. This defaults to \code{FALSE}.}
  \item{Discrep}{This optional argument allows a discrepancy statistic to
    be included. For more information on discrepancy statistics, see
    \code{summary.demonoid.ppc}.}
  \item{d}{This is an optional integer to be used with the
    \code{Discrep} argument above, and it defaults to \code{d=0}. For
    more information on discrepancy, see
    \code{\link{summary.demonoid.ppc}}.}
}
\details{Variable importance is defined here as the impact of each
  variable (predictor, or column vector) in design matrix
  \eqn{\textbf{X}}{X} on \eqn{\textbf{y}^{rep}}{y^rep} (or
  \eqn{\textbf{Y}^{rep}}{Y^rep}), when the variable is removed.

  First, the full model is predicted with the
  \code{\link{predict.demonoid}} or \code{\link{predict.laplace}}
  function, and summarized with the \code{\link{summary.demonoid.ppc}}
  or \code{\link{summary.laplace.ppc}} function, respectively. The
  results are stored in the first row of the output. Each successive
  row in the output corresponds to the application of \code{predict} and
  \code{summary} functions, but with each variable in design matrix
  \eqn{\textbf{X}}{X} being set to zero and effectively removed. The
  results show the impact of sequentially removing each predictor.

  The criteria for variable importance may differ from model to
  model. As a default, the L-criterion is recommended for models with
  dependent variables that are not categorical, provided that
  sample size is small enough that it does not result in \code{Inf}, and
  \code{Discrep="min(yhat[i,] < min(y)"} is recommended for models with
  categorical dependent variables.

  With the L-criterion, variable importance has a positive relationship, 
  such that larger values indicate a more important variable, because
  removing that variable resulted in a worse fit to the data, for
  example. Also, \code{S.L} may be helpful in interpreting the extent of
  a change in fit. Variable importance has a negative relationship with
  Mean.Lift.

  The \code{Importance} function may indicate that a model fits the data
  better with a variable removed. In which case, the user may
  choose to leave the variable in the model (perhaps the model is
  misspecified without the variable), investigate the relationship and
  re-specify the relationship between the independent and dependent
  variable(s), or remove the variable and update the model again.

  In contrast to variable importance, the \code{\link{PosteriorChecks}}
  function calculates parameter importance, which is the probability
  that each parameter's marginal posterior distribution is greater than
  zero, where an important parameter does not include zero in its
  probability interval. Parameter importance and variable importance
  may disagree, and should both be studied.

  The \code{Importance} function tends to indicate that a model fits the
  data better when variables are removed that have parameters with
  marginal posterior distributions that include 0 in the 95\%
  probability interval (variables associated with lower parameter
  importance).

  Oftentimes, in complicated models, it is difficult to assess variable
  importance by examining the marginal posterior distribution of the
  associated parameter(s). Consider polynomial regression, in which each
  variable may have multiple parameters.

  The information provided by the \code{Importance} function may be used
  for model revision, or reporting the relative importance of variables.
}
\value{
  \code{Importance} returns a matrix with a number of rows equal to the
  number of columns in design matrix \eqn{\textbf{X}}{X} + 1 (including
  the full model), and 3 columns, which are Concordance (or Mean.Lift if
  categorical), Discrep, and L-criterion. Each row represents a model
  with a predictor in \eqn{\textbf{X}}{X} removed (except for the first
  row, which is the full model), and the resulting posterior predictive
  checks.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{PosteriorChecks}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.laplace}},
  \code{\link{summary.demonoid.ppc}},
  \code{\link{summary.laplace.ppc}}.
}
\examples{
#First, update the model with the LaplacesDemon function, such as
#the example with linear regression, creating an object called Fit.
#Then
#Importance(Fit, Model, MyData, Discrep="Chi-Square")
}
\keyword{posterior predictive check}