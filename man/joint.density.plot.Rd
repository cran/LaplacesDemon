\name{joint.density.plot}
\alias{joint.density.plot}
\title{Joint Density Plot}
\description{
This function plots the joint kernel density from samples of two
marginal posterior distributions.
}
\usage{
joint.density.plot(x, y, Title=NULL, contour=TRUE, color=FALSE)
}
\arguments{
  \item{x,y}{These are vectors consisting of samples from two marginal
    posterior distributions, such as those output by
    \code{\link{LaplacesDemon}} in components \code{Posterior1} (all
    samples) or \code{Posterior2} (stationary samples).}
  \item{Title}{This is the title of the joint posterior density plot.}
  \item{contour}{This logical argument indicates whether or not contour
    lines will be added to the plot. \code{contour} defaults to
    \code{TRUE}.}
  \item{color}{This logical argument indicates whether or not color will
    be added to the plot. \code{color} defaults to \code{FALSE}.}
}
\details{
  This function produces either a bivariate scatterplot that may have
  kernel density contour lines added, or a bivariate plot with kernel 
  density-influenced colors, which may also have kernel density contour
  lines added. A joint density plot may be more informative than two
  univariate density plots.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
X <- rmvn(1000, runif(2), diag(2))
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=TRUE, color=FALSE)
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=FALSE, color=TRUE)
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=TRUE, color=TRUE)
}
\keyword{plot}