\name{dist.Wishart}
\alias{dwishart}
\alias{rwishart}
\title{The Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the Wishart distribution.
}
\usage{
   dwishart(Omega, nu, S, log=FALSE)
   rwishart(nu, S)
}
\arguments{
    \item{Omega}{This is the symmetric, positive-definite \eqn{p \times
	p}{p x p} matrix \eqn{\Omega}{Omega}.}
    \item{nu}{This is the scalar degrees of freedom \eqn{\nu}{nu}.}
    \item{S}{This is the symmetric, inverse \eqn{p \times p}{p x p} scale
        matrix \eqn{\textbf{S}}{S}.}
    \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\Omega) =
      |\textbf{S}^{-1}|^{\nu/2}|\Omega|^{(\nu-p-1)/2} \exp(-(1/2)tr[\textbf{S}^{-1}\Omega])}{p(Omega) = |S^(-1)|^(nu/2)|Omega|^[(nu-p-1)/2] exp{-(1/2)tr[S^(-1)Omega]}}
    \item Inventor: John Wishart (1928)
    \item Notation 1: \eqn{\Omega \sim \mathcal{W}_{\nu}(\textbf{S})}{Omega ~ W[nu](S)}
    \item Notation 2: \eqn{p(\Omega) = \mathcal{W}_{\nu}(\Omega | \textbf{S})}{p(Omega) = W[nu](Omega | S)}
    \item Parameter 1: degrees of freedom \eqn{\nu}{nu}
    \item Parameter 2: symmetric, non-singular \eqn{p \times p}{p x p} inverse scale matrix \eqn{\textbf{S}}{S}
    \item Mean: \eqn{E(\Omega) = \nu \textbf{S}}{E(Omega) = nuS}
    \item Variance: \eqn{var(\Omega) = \nu (\textbf{S}^2_{i,j} + \textbf{S}_{i,i} \textbf{S}_{j,j})}{var(Omega) = nu(S[i,j]^2 + S[i,i]S[j,j])}
    \item Mode: \eqn{mode(\Omega) = (\nu - p - 1)
      \textbf{S}}{mode(Omega) = (nu-p-1)S}, for \eqn{\nu \ge p + 1}{nu >= p + 1}
}

The Wishart distribution is a generalization to multiple dimensions of
the chi-square distribution, or, in the case of non-integer degrees of
freedom, of the gamma distribution. The Wishart is the conjugate prior
distribution for the inverse covariance matrix (the precision matrix)
\eqn{\Omega}{Omega}, the inverse of which (covariance matrix
\eqn{\Sigma}{Sigma}) is used in a multivariate normal distribution.
The integral is finite if the scalar degrees of freedom parameter,
\eqn{\nu}{nu}, is greater than or equal to the dimension \eqn{p} of
the symmetric, inverse scale matrix \eqn{\textbf{S}}{S}. The degrees
of freedom, \eqn{\nu}{nu}, is equivalent to specifying a prior sample
size, indicating the confidence in \eqn{\textbf{S}}{S}, where
\eqn{\textbf{S}}{S} is a prior guess at the order of covariance matrix
\eqn{\Sigma}{Sigma}. An uninformative prior distribution is obtained
as \eqn{\nu \rightarrow 0}{nu -> 0}.

These functions are similar to those in the \code{MCMCpack} package.
}
\value{
  \code{dwishart} gives the density and
  \code{rwishart} generates random deviates.
}
\examples{
x <- dwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{distribution}
