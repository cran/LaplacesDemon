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
    \item{Sigma}{This is the symmetric, positive-definite
         \eqn{p \times p}{p x p} matrix \eqn{\Sigma}{Sigma}.}
    \item{nu}{This is the scalar degrees of freedom.}
    \item{R}{This is the \eqn{p \times p}{p x p} symmetric scale
         matrix \eqn{\textbf{R}}{R}.}
    \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \deqn{p(\Sigma) =
      \frac{|\textbf{R}|^{\nu/2}|\Sigma|^{-(\nu + p + 1)} \exp(-tr(\textbf{R}\Sigma^{-1})/2)}{2^{\nu p / 2} \Gamma_p(\nu/2)}}{p(Sigma) = [|R|^(nu/2)|Sigma|^[-(nu+p+1)/2]exp(-tr[R*Sigma^(-1)]/2)]/[2^(nu*p/2)Gamma[p](nu/2)]}
    \item Inventor: John Wishart (1928)
    \item Notation 1: \eqn{\Sigma \sim \mathcal{W}^{-1}_{\nu}(\textbf{R})}{Sigma ~ W^-1[\code{nu}](\code{R})}
    \item Notation 2: \eqn{p(\Sigma) = \mathcal{W}^{-1}_{\nu}(\Sigma | \textbf{R})}{p(Sigma) = W^-1[\code{nu}](\code{Sigma} | \code{R})}
    \item Parameter 1: degrees of freedom \eqn{\nu}{nu}
    \item Parameter 2: symmetric, non-singular \eqn{p \times p}{p x p} scale matrix \eqn{\textbf{R}}{R}
    \item Mean: \eqn{E(\Sigma) = \frac{\textbf{R}}{\nu - p - 1}}{E(Sigma) = R / (nu - p - 1)}
    \item Variance: 
    \item Mode: \eqn{mode(\Sigma) = \frac{\textbf{R}}{\nu + p + 1}}{mode(Sigma) = R / (nu + p + 1)}
}

The inverse Wishart distribution is a probability distribution defined on
real-valued positive-definite matrices, and is used as the conjugate
prior for the covariance matrix, \eqn{\Sigma}{Sigma}, of a multivariate
normal distribution. The inverse-Wishart density is always finite, and the
integral is always finite. A degenerate form occurs when \eqn{\nu < p}{nu < p}.

These functions are similar to those in the \code{MCMCpack} package.
}
\value{
  \code{dinvwishart} gives the density and
  \code{rinvwishart} generates random deviates.
}
\examples{
x <- dinvwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rinvwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}
