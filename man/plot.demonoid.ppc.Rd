\name{plot.demonoid.ppc}
\alias{plot.demonoid.ppc}
\title{Plots of Posterior Predictive Checks}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{demonoid.ppc}. A variety of plots is provided.
}
\usage{\method{plot}{demonoid.ppc}(x, Style=NULL, Data=NULL, Rows=NULL,
     PDF=FALSE, \dots)}
\arguments{
     \item{x}{
     This required argument is an object of class \code{demonoid.ppc}.
     } 
     \item{Style}{
     This optional argument specifies one of several styles of plots,
     and defaults to \code{NULL} (which is the same as
     \code{"Density"}). Styles of plots are indicated in quotes. Optional
     styles include \code{"Covariates"}, \code{"Covariates, Categorical DV"},
     \code{"Density"}, \code{"Fitted"}, \code{"Fitted, Multivariate, C"},
     \code{"Fitted, Multivariate, R"}, \code{"Predictive Quantiles"},
     \code{"Residuals"}, \code{"Residuals, Multivariate, C"},
     \code{"Residuals, Multivariate, R"}, \code{"Space-Time by Space"},
     \code{"Space-Time by Time"}, \code{"Spatial"},
     \code{"Spatial Uncertainty"}, \code{"Time-Series"},
     \code{"Time-Series, Multivariate, C"}, and
     \code{"Time-Series, Multivariate, R"}. Details are given below.
     }
     \item{Data}{
     This optional argument accepts the data set used when updating the
     model. Data is required only with certain plot styles, including
     \code{"Covariates"}, \code{"Covariates, Categorical DV"},
     \code{"Fitted, Multivariate, C"}, \code{"Fitted, Multivariate, R"},
     \code{"Residuals, Multivariate, C"}, \code{"Residuals, Multivariate, R"},
     \code{"Space-Time by Space"}, \code{"Space-Time by Time"},
     \code{"Spatial"}, \code{"Spatial Uncertainty"},
     \code{"Time-Series, Multivariate, C"}, and
     \code{"Time-Series, Multivariate, R"}.
     }
     \item{Rows}{
     This optional argument is for a vector of row numbers that 
     specify the records associated by row in the object of class
     \code{demonoid.ppc}. Only these rows are plotted. The default is to
     plot all rows. Some plots do not allow rows to be specified.
     }
     \item{PDF}{
     This logical argument indicates whether or not the user wants
     Laplace's Demon to save the plots as a .pdf file.
     }
     \item{\dots}{Additional arguments are unused.}
}
\details{
  This function can be used to produce a variety of posterior predictive
  plots, and the style of plot is selected with the \code{Style}
  argument. Below are some notes on the styles of plots.

  \code{Covariates} requires \code{Data} to be specified, and also
  requires that the covariates are named \code{X} or \code{x}. A plot
  is produced for each covariate column vector against yhat, and is
  appropriate when y is not categorical.

  \code{Covariates, Categorical DV} requires \code{Data} to be
  specified, and also requires that the covariates are named \code{X} or
  \code{x}. A plot is produced for each covariate column vector against
  yhat, and is appropriate when y is categorical.
  
  \code{Density} plots show the kernel density of the posterior
  predictive distribution for each selected row of y (all are selected
  by default). A vertical red line indicates the position of the
  observed y along the x-axis. When the vertical red line is close to
  the middle of a normal posterior predictive distribution, then there
  is little discrepancy between y and the posterior predictive
  distribution. When the vertical red line is in the tail of the
  distribution, or outside of the kernel density altogether, then
  there is a large discrepancy between y and the posterior predictive
  distribution. Large discrepancies may be considered outliers, and
  moreover suggest that an improvement in model fit should be considered.

  \code{Fitted} plots compare y with its replicate and provide loess
  smoothing. Otherwise, y is also compared with the quantiles of the
  probability interval of its replicate, and smoothing is also
  applied. This plot is appropriate when y is univariate and not
  categorical.

  \code{Fitted, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each column-wise
  vector of y in Y with its replicate and provide loess smoothing.
  Otherwise, y is also compared with the quantiles of the probability
  interval of its replicate, and smoothing is also applied. This plot
  is appropriate when Y is multivariate, not categorical, and desired to
  be seen column-wise.
  
  \code{Fitted, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each row-wise
  vector of y in Y with its replicate and provide loess smoothing.
  Otherwise, y is also compared with the quantiles of the probability
  interval of its replicate, and smoothing is also applied. This plot
  is appropriate when Y is multivariate, not categorical, and desired to
  be seen row-wise.
  
  \code{Predictive Quantiles} plots compare y with the predictive
  quantile (PQ) of its replicate. This may be useful in looking for
  patterns with outliers. Instances outside of the gray lines are
  considered outliers. 

  \code{Residuals} plots compare y with its residuals. The probability
  interval is plotted as a line in gray. This plot is appropriate when y
  is univariate.

  \code{Residuals, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are plots of each column-wise
  vector of residuals. The probability interval is plotted as a line in
  gray. This plot is appropriate when Y is multivariate, not
  categorical, and the residuals are desired to be seen column-wise.

  \code{Residuals, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are plots of each row-wise
  vector of residuals. The probability interval is plotted as a line in
  gray. This plot is appropriate when Y is multivariate, not
  categorical, and the residuals are desired to be seen row-wise.
  
  \code{Space-Time by Space} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude}, \code{longitude}, \code{S}, and
  \code{T}. These space-time plots compare the S x T matrix Y with the S
  x T matrix Yrep, producing one time-series plot per point s in space,
  for a total of S plots. Therefore, these are time-series plots for
  each point s in space across T time-periods. See \code{Time-Series}
  plots below.

  \code{Space-Time by Time} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude}, \code{longitude}, \code{S}, and
  \code{T}. These space-time plots compare the S x T matrix Y with the S
  x T matrix Yrep, producing one spatial plot per time-period, and T
  plots will be produced. See \code{Spatial} plots below.

  \code{Spatial} requires \code{Data} to be specified, and also requires
  that the following variables exist in the data set with exactly these
  names: \code{latitude} and \code{longitude}. This spatial plot shows
  yrep plotted according to its coordinates, and is color-coded so that
  higher values of yrep become more red, and lower values become more
  yellow.

  \code{Spatial Uncertainty} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude} and \code{longitude}. This
  spatial plot shows the probability interval of yrep plotted according
  to its coordinates, and is color-coded so that wider probability
  intervals become more red, and lower values become more yellow.

  \code{Time-Series} plots compare y with its replicate, including the
  median and probability interval quantiles. This plot is appropriate
  when y is univariate and ordered by time.

  \code{Time-Series, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each column-wise
  time-series in Y with its replicate, including the median and
  probability interval quantiles. This plot is appropriate when y is
  multivariate and each time-series is indexed by column in Y.

  \code{Time-Series, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each row-wise
  time-series in Y with its replicate, including the median and
  probability interval quantiles. This plot is appropriate when y is
  multivariate and each time-series is indexed by row in Y, such as is
  typically true in panel models.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{predict.demonoid}}}
\examples{### See the LaplacesDemon function for an example.}
\keyword{plot}