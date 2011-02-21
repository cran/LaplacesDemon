\name{dist.Dirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\title{The Dirichlet Distribution}
\description{
  This is the density function and random generation from the Dirichlet
  distribution.
}
\usage{
ddirichlet(x, alpha, log=FALSE)
rdirichlet(n, alpha)
}
\arguments{
  \item{x}{This is a vector containing a single deviate or matrix
    containing one random deviate per row.}
  \item{n}{This is the number of random vectors to generate.}
  \item{alpha}{This is a vector of shape parameters, or matrix of shape
    parameters corresponding to the number of draws.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: p(theta) = (gamma(alpha[1] + ... + alpha[k]) / 
          (gamma(alpha[1]) ... gamma(alpha[k]))) * 
          theta[1]^(alpha[1] - 1) ... theta[k]^(alpha[k] - 1)
          theta[1],...,theta[k] >= 0, 
          (the sum of j=1 to k of) theta[j] = 1
    \item Inventor: Johann Peter Gustav Lejeune Dirichlet (1805-1859)
    \item Notation 1: theta ~ Dirichlet(alpha[1],..., alpha[k])
    \item Notation 2: p(theta) = Dirichlet(theta | alpha[1],...,alpha[k])
    \item Notation 3: theta ~ Dir(alpha[1],..., alpha[k])
    \item Notation 4: p(theta) = Dir(theta | alpha[1],...,alpha[k])
    \item Parameter: 'prior sample sizes' alpha[j] > 0, alpha[0] = (the sum from j=1 to k of) alpha[j]
    \item Mean: E(theta[j]) = alpha[j] / alpha[0]
    \item Variance: var(theta[j]) = (alpha[j] * (alpha[0] - alpha[j])) / 
    (alpha[0]^2 * (alpha[0] + 1))
    \item Covariance: cov(theta[i], theta[j]) = - ((alpha[i]*alpha[j]) / 
          (alpha[0]^2 * (alpha[0] + 1)))
    \item Mode: mode(theta[j]) = (alpha[j] - 1) / (alpha[0] - k)
}

The Dirichlet is the conjugate prior distribution for the parameters of
the categorical and multinomial distributions. The Dirichlet is a 
multivariate generalization of the beta distribution. Its probability
density function returns the belief that the probabilities of k rival
events are theta[j] given that each event has been observed alpha[j] - 1
times.

These functions are similar to those in the \code{MCMCpack} package.
}
\value{
  \code{ddirichlet} gives the density. \code{rdirichlet} returns a
  matrix with \code{n} rows, each containing a single Dirichlet random
  deviate.
}
\examples{
x <- ddirichlet(c(.1,.3,.6), c(1,1,1))
x <- rdirichlet(10, c(1,1,1))
}
\keyword{distribution}
