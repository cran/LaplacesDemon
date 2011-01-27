\name{invwishart}
\alias{dinvwishart}
\alias{rinvwishart}
\title{The Inverse Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the inverse Wishart distribution.
}
\usage{
   dinvwishart(W, v, S)
   rinvwishart(v, S)
}
\arguments{
    \item{W}{This is positive-definite matrix W \eqn{(p \times p)}{(p x p)}.}
    \item{v}{This is the scalar degrees of freedom, and is more
      correctly called \code{nu} (which looks like a \code{v}).}
    \item{S}{This is the scale matrix \eqn{(p \times p)}{(p x p)}.}
}
\value{
  \code{dinvwishart} evaluates the density at positive definite matrix \code{W}.
  \code{rinvwishart} generates one random draw from the distribution.
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: p(W) = [|S|^(v/2)det(W)^[-(v+k+1)/2]exp(-Tr[SW^(-1)]/2)]/[2^(vk/2)Gamma[k](v/2)]
    \item Inventor: John Wishart (1928)
    \item Notation 1: IW ~ Inv-Wishart[\code{v}](\code{S})
    \item Notation 2: p(IW) = Inv-Wishart[\code{v}](\code{W} | \code{S})
    \item Parameter 1: degrees of freedom \code{v} (or nu)
    \item Parameter 2: symmetric, positive-definite k x k scale matrix \code{S}
    \item Mean: E(W) = S / (v - k - 1)
    \item Variance: 
    \item Mode: S / (v + k + 1)
}

The inverse Wishart distribution is a probability distribution defined on
real-valued positive-definite matrices, and is used as the conjugate
prior for the covariance matrix of a multivariate normal distribution.

These functions are similar to those in the \code{MCMCpack} package.
}
\examples{
x <- dinvwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rinvwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}



