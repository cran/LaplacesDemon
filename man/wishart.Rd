\name{wishart}
\alias{dwishart}
\alias{rwishart}
\title{The Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the Wishart distribution.
}
\usage{
   dwishart(W, v, S)
   rwishart(v, S)
}
\arguments{
    \item{W}{This is positive-definite matrix W \eqn{(p \times p)}{(p x p)}.}
    \item{v}{This is the scalar degrees of freedom, and is more
      correctly called \code{nu} (which looks like a \code{v}).}
    \item{S}{This is the inverse scale matrix \eqn{(p \times p)}{(p x p)}.}
}
\value{
  \code{dwishart} evaluates the density at positive definite matrix \code{W}.
  \code{rwishart} generates one random draw from the distribution.
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: p(W) = ({|W|^[(n-k-1)/2]}/[2^(nk/2)det(V)^(n/2)Gamma[k](n/2)])exp[0.5Tr(V^-1W)]
    \item Inventor: John Wishart (1928)
    \item Notation 1: W ~ Wishart[\code{v}](\code{S})
    \item Notation 2: p(W) = Wishart[\code{v}](\code{W} | \code{S})
    \item Parameter 1: degrees of freedom \code{v} (or nu)
    \item Parameter 2: symmetric, positive-definite k x k scale matrix \code{S}
    \item Mean: E(W) = vS
    \item Variance: n(v[i,j]^2 + v[i,i]v[j,j])
    \item Mode: (n-k-1)V, for n >= k + 1
}

The Wishart distribution is a generalization to multiple dimensions of
the chi-square distribution, or, in the case of non-integer degrees of
freedom, of the gamma distribution. The Wishart is the conjugate prior
distribution for the inverse covariance matrix (the precision matrix) in
a multivariate normal distribution.  The integral is finite if the
degrees of freedom parameter, \code{v}, is greater than or equal to the
dimension, \code{k}. An uninformative distribution is obtained as
\code{v -> 0}.

These functions are similar to those in the \code{MCMCpack} package.
}
\examples{
x <- dwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}
