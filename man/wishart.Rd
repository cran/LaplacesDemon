\name{wishart}
\alias{dwishart}
\alias{rwishart}
\title{The Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the Wishart distribution.
}
\usage{
   dwishart(Omega, nu, S, log = FALSE)
   rwishart(nu, S)
}
\arguments{
    \item{Omega}{This is symmetric, positive-definite matrix Omega \eqn{(p \times p)}{(p x p)}.}
    \item{nu}{This is the scalar degrees of freedom}
    \item{S}{This is the symmetric inverse scale matrix \eqn{(p \times p)}{(p x
	p)}.}
    \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\value{
  \code{dwishart} evaluates the density at positive definite matrix \code{Omega}.
  \code{rwishart} generates one random draw from the distribution.
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: p(Omega) = |S^(-1)|^(nu/2)|Omega|^[(nu-p-1)/2] exp{-(1/2)tr[S^(-1)Omega]}
    \item Inventor: John Wishart (1928)
    \item Notation 1: Omega ~ Wishart[\code{nu}](\code{S})
    \item Notation 2: p(Omega) = Wishart[\code{nu}](\code{Omega} | \code{S})
    \item Parameter 1: degrees of freedom \code{nu}
    \item Parameter 2: symmetric, non-singular p x p inverse scale matrix \code{S}
    \item Mean: E(Omega) = nuS
    \item Variance: nu(S[i,j]^2 + S[i,i]S[j,j])
    \item Mode: (nu-p-1)S, for nu >= p + 1
}

The Wishart distribution is a generalization to multiple dimensions of
the chi-square distribution, or, in the case of non-integer degrees of
freedom, of the gamma distribution. The Wishart is the conjugate prior
distribution for the inverse covariance matrix (the precision matrix)
\code{Omega}, the inverse of which (covariance matrix \code{Sigma}) is
used in a multivariate normal distribution.  The integral is finite if the
scalar degrees of freedom parameter, \code{nu}, is greater than or equal
to the dimension \code{p} of the symmetric, inverse scale matrix
\code{S}. The degrees of freedom, \code{nu}, is equivalent to specifying
a prior sample size, indicating the confidence in \code{S}, where
\code{S} is a prior guess at the order of covariance matrix
\code{Sigma}. An uninformative distribution is obtained as \code{nu -> 0}.

These functions are similar to those in the \code{MCMCpack} package.
}
\examples{
x <- dwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}
