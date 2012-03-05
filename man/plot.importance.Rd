\name{plot.importance}
\alias{plot.importance}
\title{Plot Variable Importance}
\description{
  This may be used to plot variable importance with the L-criterion
  regarding an object of class \code{importance}.
}
\usage{\method{plot}{importance}(x, \dots)}
\arguments{
  \item{x}{
    This required argument is an object of class \code{importance}.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The x-axis is the L-criterion of the \code{\link{Importance}}
  function, and variables are on the y-axis. A more important variable
  is associated with a dot that is plotted farther to the right. For
  more information on variable importance, see the
  \code{\link{Importance}} function.
}
\references{
  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{Importance}}}
\keyword{plot}