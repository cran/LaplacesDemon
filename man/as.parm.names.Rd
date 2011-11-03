\name{as.parm.names}
\alias{as.parm.names}
\title{Parameter Names}
\description{
This function creates a vector of parameter names from a list of
parameters, and the list may contain any combination of scalars,
vectors, matrices, upper-triangular matrices, and arrays.
}
\usage{
as.parm.names(x, uppertri=NULL)
}
\arguments{
     \item{x}{This required argument is a list of named parameters. The
       list may contain scalars, vectors, matrices, and arrays.}
     \item{uppertri}{This optional argument must be a vector with a
       length equal to the number of named parameters. Each element in
       \code{uppertri} must be either a 0 or 1, where a 1 indicates that
       an upper triangular matrix will be used for the associated
       element in the vector of named parameters. Each element of
       \code{uppertri} is associated with a named parameter. The
       \code{uppertri} argument does not function with arrays.}
}
\details{
Each \code{model} function for \code{\link{LaplaceApproximation}} or
\code{\link{LaplacesDemon}} requires a vector of parameters (specified
at first as \code{Initial.Values}) and a list of data. One component in
the list of data must be named parm.names. Each element of parm.names
is a name associated with the corresponding parameter in
\code{Initial.Values}.

The parm.names vector is easy to program explicitly for a simple model,
but can require considerably more programming effort for more
complicated models. The \code{as.parm.names} function is a utility
function designed to minimize programming by the user.

For example, a simple model may only require \code{parm.names <-
  c("alpha", "beta[1]", "beta[2]", "log.sigma")}. A more complicated model
  may contain hundreds of parameters that are combinations of scalars,
  vectors, matrices, upper-triangular matrices, and arrays, and is the
  reason for the \code{as.parm.names} function.

In the case of an upper-triangular matrix, simply pass the full matrix
to \code{as.parm.names} and indicate that only the upper-triangular will
be used via the \code{uppertri} argument.

Numerous examples may be found in the accompanying ``Examples'' vignette.
}
\value{
 This function returns a vector of parameter names.
}
\references{
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\examples{
N <- 100
J <- 5
y <- rnorm(N,0,1)
X <- matrix(runif(N*J,-2,2),N,J)
S <- diag(J)
T <- diag(2)
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(log.sigma=0, beta=rep(0,J), S=diag(J),
     T=diag(2)), uppertri=c(0,0,0,1))
MyData <- list(J=J, N=N, S=S, T=T, X=X, mon.names=mon.names,
     parm.names=parm.names, y=y)
MyData
}
\keyword{utility}