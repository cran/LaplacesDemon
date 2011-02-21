\name{dist.Inverse.Wishart}
\alias{dinvwishart}
\alias{rinvwishart}
\title{The Inverse Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the inverse Wishart distribution.
}
\usage{
   dinvwishart(Sigma, nu, R, log=FALSE)
   rinvwishart(nu, R)
}
\arguments{
    \item{Sigma}{This is the symmetric, positive-definite matrix Sigma \eqn{(p \times p)}{(p x p)}.}
    \item{nu}{This is the scalar degrees of freedom.}
    \item{R}{This is the symmetric scale matrix \eqn{(p \times p)}{(p x p)}.}
    \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\value{
  \code{dinvwishart} evaluates the density at symmetric, positive-definite matrix \code{Sigma}.
  \code{rinvwishart} generates one random draw from the distribution.
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: p(Sigma) = [|R|^(nu/2)|Sigma|^[-(nu+p+1)/2]exp(-tr[R*Sigma^(-1)]/2)]/[2^(nu*p/2)Gamma[p](nu/2)]
    \item Inventor: John Wishart (1928)
    \item Notation 1: Sigma ~ Inv-Wishart[\code{nu}](\code{R})
    \item Notation 2: p(Sigma) = Inv-Wishart[\code{nu}](\code{Sigma} | \code{R})
    \item Parameter 1: degrees of freedom \code{nu}
    \item Parameter 2: symmetric, non-singular p x p scale matrix \code{R}
    \item Mean: E(Sigma) = R / (nu - p - 1)
    \item Variance: 
    \item Mode: R / (nu + p + 1)
}

The inverse Wishart distribution is a probability distribution defined on
real-valued positive-definite matrices, and is used as the conjugate
prior for the covariance matrix, \code{Sigma}, of a multivariate normal
distribution. The inverse-Wishart density is always finite, and the
integral is always finite. A degenerate form occurs when \code{nu < p}.

These functions are similar to those in the \code{MCMCpack} package.
}
\examples{
x <- dinvwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rinvwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}
