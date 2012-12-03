\name{is.data}
\alias{is.data}
\title{is.data}
\description{
  This function provides a logical test of whether or not a given list
  of data meets minimum criteria to be considered data for
  \code{\link{LaplaceApproximation}} or \code{\link{LaplacesDemon}}.
}
\usage{
is.data(Data)
}
\arguments{
  \item{Data}{This is a list of data. For more information, see the
    \code{\link{LaplacesDemon}} function.}
}
\details{
  This function is useful for testing whether or not a list of data
  meets minimum criteria to be considered data in this package. The
  minimum requirements are that \code{Data} is a list, and it contains
  \code{mon.names} and \code{parm.names}.

  This function is not extensive. For example, it does not match the
  length of \code{parm.names} with the length of \code{Initial.Values},
  or compare the length of \code{mon.names} to the number of monitored
  variables output from the \code{Model} specification
  function. Additional checks are conducted in
  \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}}.
}
\value{
  The \code{is.data} function returns a logical value. It returns
  \code{TRUE} if \code{Data} meets minimum requirements to be considered
  data in this package, and \code{FALSE} otherwise.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplaceApproximation}} and
  \code{\link{LaplacesDemon}}.
}
\keyword{Utility}