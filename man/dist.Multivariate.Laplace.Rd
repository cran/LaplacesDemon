\name{dist.Multivariate.Laplace}
\alias{dmvl}
\alias{rmvl}
\title{Multivariate Laplace Distribution}
\description{
  These functions provide the density and random number generation
  for the multivariate Laplace distribution.
}
\usage{
dmvl(x, mu, Sigma, log=FALSE)
rmvl(n, mu, Sigma)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{Sigma}{This is the \eqn{k \times k}{k x k} covariance matrix
       \eqn{\Sigma}{Sigma}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) = \frac{k\Gamma(k/2)}{\pi^{k/2}
      \sqrt{|\Sigma|} \Gamma(1 + k) 2^{1 + k}}
    \exp(-\frac{1}{2}(\theta-\mu)^T \Sigma
    (\theta-\mu))^{1/2}}{p(theta) = ((k*Gamma(k/2)) / (pi^(k/2) * 
    sqrt(|Sigma|) * Gamma(1 + k) * 2^(1 + k))) *
          exp(-(1/2)*(theta-mu)^T Sigma(theta-mu))^(1/2)}
  \item Inventor: Gomez, Gomez-Villegas, and Marin (1998)
  \item Notation 1: \eqn{\theta \sim \mathcal{MVL}(\mu, \Sigma)}{theta ~
    MVL(mu, Sigma)}
  \item Notation 2: \eqn{\theta \sim \mathcal{L}_k(\mu, \Sigma)}{theta ~
    L[k](mu, Sigma)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MVL}(\theta | \mu,
    \Sigma)}{p(theta) = MVL(theta | mu, Sigma)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{L}_k(\theta | \mu,
    \Sigma)}{p(theta) = L[k](theta | mu, Sigma)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k}
  covariance matrix \eqn{\Sigma}{Sigma}
  \item Mean: \eqn{E(\theta) = }{E(theta) = }
  \item Variance: \eqn{var(\theta) =}{var(theta) = }
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The multivariate Laplace distribution is a multidimensional extension of
the one-dimensional or univariate symmetric Laplace distribution. There
are multiple formulas for the multivariate Laplace distribution. This
multivariate Laplace distribution is the one that is expressed as a
reduced form of the multivariate power exponential distribution, when
its kurtosis parameter is \eqn{\kappa = 0.5}{kappa = 0.5}.

The multivariate Laplace distribution is an attractive alternative to
the multivariate normal distribution due to its wider tails, and remains
a two-parameter distribution, unlike the three-parameter multivariate t
distribution, which is often used as a robust alternative to the
multivariate normal distribution.
}
\value{
  \code{dmvl} gives the density, and
  \code{rmvl} generates random deviates.
}
\references{
  Gomez, E., Gomez-Villegas, M.A., and Marin, J.M. (1998). "A
  Multivariate Generalization of the Power Exponential Family of
  Distributions". Communications in Statistics-Theory and Methods, 27(3),
  p. 589--600.
}
\author{Statisticat, LLC. \email{statisticat@gmail.com}}
\seealso{
  \code{\link{dlaplace}},
  \code{\link{dmvn}},
  \code{\link{dmvnp}},
  \code{\link{dmvpe}},
  \code{\link{dmvt}},
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
n <- 100
k <- 3
x <- matrix(runif(n*k),n,k)
mu <- matrix(runif(n*k),n,k)
Sigma <- diag(k)
dmvl(x, mu, Sigma)
}
\keyword{Distribution}