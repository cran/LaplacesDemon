\name{bernoulli}
\alias{dbern}
\alias{pbern}
\alias{qbern}
\alias{rbern}
\title{The Bernoulli Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Bernoulli distribution with parameter \code{prob}.
}
\usage{
dbern(x, prob, log = FALSE)
pbern(q, prob, lower.tail = TRUE, log.p = FALSE)
qbern(p, prob, lower.tail = TRUE, log.p = FALSE)
rbern(n, prob)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{prob}{probability of success on each trial.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dbern} gives the density, \code{pbern} gives the distribution
  function, \code{qbern} gives the quantile function and \code{rbern}
  generates random deviates.
}
\details{
  The Bernoulli distribution with \code{prob} \eqn{= p} has density
  \deqn{p(x) = {p}^{x} {(1-p)}^{1-x}}{
    p(x) = p^x (1-p)^(1-x)}
  for \eqn{x = 0 or 1}.

  If an element of \code{x} is not \code{0} or \code{1}, the result of \code{dbern}
  is zero, without a warning.
  \eqn{p(x)} is computed using Loader's algorithm, see the reference below.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\seealso{
  \code{\link{dbinom}} for the binomial (Bernoulli is a special case
  of the binomial), and \code{\link{dpois}} for the Poisson distribution.
}
\examples{
# Compute P(X=1) for X Bernoulli(0.7)
dbern(1, 0.7)
}
\keyword{distribution}
