\name{CSF}
\alias{CSF}
\title{Cumulative Sample Function}
\description{
  The Cumulative Sample Function (CSF) is a visual MCMC diagnostic in
  which the user may select a measure (such as a variable, summary
  statistic, or other diagnostic), and observe a plot of how the measure
  changes over cumulative posterior samples from MCMC, such as the
  output of \code{\link{LaplacesDemon}}. As a generalized visual
  diagnostic across a variety of measures, this function is believed to
  be a novel contribution to MCMC diagnostics, introduced in the
  LaplacesDemon package. This may be considered to be a generalized
  extension of the \code{cumuplot} in the coda package.
}
\usage{
CSF(x, name, method="Quantiles", quantiles=c(0.025,0.500,0.975), output=FALSE)
}
\arguments{
  \item{x}{This is a vector of posterior samples from MCMC.}
  \item{name}{This is an optional name for vector \code{x}, and is input
    as a quoted string, such as \code{name="theta"}.}
  \item{method}{This is a measure that will be observed over the course
    of cumulative samples of \code{x}. It defaults to
    \code{method="Quantiles"}, and optional methods include: \code{"ESS"},
    \code{"Geweke.Diagnostic"}, \code{"HPD"}, \code{"is.stationary"},
    \code{"Kurtosis"}, \code{"MCSE"}, \code{"MCSE.bm"},
    \code{"MCSE.sv"}, \code{"Mean"}, \code{"Mode"}, \code{"N.Modes"},
    \code{"Precision"}, \code{"Quantiles"}, and \code{"Skewness"}.}
  \item{quantiles}{This optional argument applies only when
    \code{method="Quantiles"}, in which case this vector indicates the
    probabilities that will be observed. It defaults to the median and
    95\% probability interval bounds.}
  \item{output}{Logical. If \code{output=TRUE}, then the results of the
    measure over the course of the cumulative samples will be output as
    an object, either a vector or matrix, depending on the \code{method}
  argument. The \code{output} argument defaults to \code{FALSE}.}
}
\details{
  When \code{method="ESS"}, the effective sample size (ESS) is observed
  as a function of the cumulative samples of \code{x}. For more
  information, see the \code{\link{ESS}} function.

  When \code{method="Geweke.Diagnostic"}, the Z-score output of the
  Geweke diagnostic is observed as a function of the cumulative samples
  of \code{x}. For more information, see the
  \code{\link{Geweke.Diagnostic}} function.

  When \code{method="HPD"}, the Highest Posterior Density (HPD) interval
  is observed as a function of the cumulative samples of \code{x}. For
  more information, see the \code{\link{HPD}} function.

  When \code{method="is.stationary"}, stationarity is logically
  tested and the result is observed as a function of the cumulative
  samples of \code{x}. For more information, see the
  \code{\link{is.stationary}} function.

  When \code{method="Kurtosis"}, kurtosis is observed as a function of
  the cumulative samples of \code{x}.

  When \code{method="MCSE"}, the Monte Carlo Standard Error (MCSE)
  estimated with the \code{IMPS} method is observed as a function of
  the cumulative samples of \code{x}. For more information, see the
  \code{\link{MCSE}} function.

  When \code{method="MCSE.bm"}, the Monte Carlo Standard Error (MCSE)
  estimated with the \code{batch.means} method is observed as a
  function of the cumulative samples of \code{x}. For more information,
  see the \code{\link{MCSE}} function.

  When \code{method="MCSE.sv"}, the Monte Carlo Standard Error (MCSE)
  estimated with the \code{sample.variance} method is observed as a
  function of the cumulative samples of \code{x}. For more information,
  see the \code{\link{MCSE}} function.

  When \code{method="Mean"}, the mean is observed as a function of
  the cumulative samples of \code{x}.

  When \code{method="Mode"}, the estimated mode is observed as a
  function of the cumulative samples of \code{x}. For more information,
  see the \code{\link{Mode}} function.

  When \code{method="N.Modes"}, the estimated number of modes is
  observed as a function of the cumulative samples of \code{x}. For
  more information, see the \code{\link{Modes}} function.

  When \code{method="Precision"}, the precision (inverse variance) is
  observed as a function of the cumulative samples of \code{x}.

  When \code{method="Quantiles"}, the quantiles selected with the
  \code{quantiles} argument are observed as a function of the
  cumulative samples of \code{x}.

  When \code{method="Skewness"}, skewness is observed as a function of
  the cumulative samples of \code{x}.
}
\references{
  \url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{ESS}},
  \code{\link{Geweke.Diagnostic}},
  \code{\link{HPD}},
  \code{\link{is.stationary}},
  \code{\link{LaplacesDemon}},
  \code{\link{MCSE}},
  \code{\link{Mode}},
  \code{\link{Modes}}.
}
\examples{
#Commented-out because of run-time for package builds
#library(LaplacesDemon)
#x <- rnorm(1000)
#CSF(x, method="ESS")
#CSF(x, method="Geweke.Diagnostic")
#CSF(x, method="HPD")
#CSF(x, method="is.stationary")
#CSF(x, method="Kurtosis")
#CSF(x, method="MCSE")
#CSF(x, method="MCSE.bm")
#CSF(x, method="MCSE.sv")
#CSF(x, method="Mean")
#CSF(x, method="Mode")
#CSF(x, method="N.Modes")
#CSF(x, method="Precision")
#CSF(x, method="Quantiles")
#CSF(x, method="Skewness")
}
\keyword{diagnostic, utility}
