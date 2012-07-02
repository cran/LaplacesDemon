\name{SIR}
\alias{SIR}
\title{Sampling Importance Resampling}
\description{
  The \code{SIR} function performs Sampling Importance Resampling, also
  called Sequential Importance Resampling, and uses a multivariate normal
  proposal density.
}
\usage{
SIR(Model, Data, mu, Sigma, n=1000)
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplaceApproximation}}.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplaceApproximation}}.}
  \item{mu}{This is a mean vector, \eqn{\mu}{mu}, for a multivariate
    normal distribution, and is usually the posterior modes from an
    object of class \code{laplace}, output from 
    \code{\link{LaplaceApproximation}}.}
  \item{Sigma}{This is a covariance matrix, \eqn{\Sigma}{Sigma}, for a
    multivariate normal distribution, and is usually the \code{Covar}
    component of an object of class \code{laplace}, output from
    \code{\link{LaplaceApproximation}}.}
  \item{n}{This is the number of samples to be drawn from the posterior
    distribution.}
}
\details{
  Sampling Importance Resampling (SIR) was introduced in Gordon, et
  al. (1993), and is the original particle filtering algorithm (and this
  family of algorithms is also known as Sequential Monte Carlo). A
  distribution is approximated with importance weights, which are
  approximations to the relative posterior densities of the particles,
  and the sum of the weights is one. In this terminology, each sample in
  the distribution is a ``particle''. SIR is a sequential or recursive
  form of importance sampling. As in importance sampling, the
  expectation of a function can be approximated as a weighted
  average. The optimal proposal distribution is the target distribution.

  In the \code{LaplacesDemon} package, the main use of the \code{SIR}
  function is to produce posterior samples for Laplace Approximation,
  and \code{SIR} is called behind-the-scenes by the
  \code{\link{LaplaceApproximation}} function. Laplace Approximation 
  estimates the posterior mode and the associated covariance
  matrix. Assuming normality, this output characterizes the marginal
  posterior distributions. However, it is often useful to have posterior
  samples, in which case the \code{SIR} function is used to draw
  samples. The number of samples, \code{n}, should increase with the
  number and intercorrelations of the parameters. Otherwise, multimodal
  posterior distributions may occur.

  In other contexts (for which this function is not designed), SIR is
  used with dynamic linear models (DLMs) and state-space models (SSMs)
  for state filtering.
  
  This function is similar to the \code{sir} function in the
  \code{LearnBayes} package.
}
\value{
  The \code{SIR} function returns a matrix of samples drawn from the
  posterior distribution.
}
\references{
  Gordon, N.J., Salmond, D.J., and Smith, A.F.M. (1993). "Novel Approach
  to Nonlinear/Non-Gaussian Bayesian State Estimation". IEEE Proceedings
  F on Radar and Signal Processing, 140 (2), p. 107--113.

  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{dmvn}} and 
  \code{\link{LaplaceApproximation}}.
}
\keyword{Sampling Importance Resampling}
\keyword{Sequential Importance Resampling}
\keyword{Utility}
