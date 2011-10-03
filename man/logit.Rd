\name{logit}
\alias{invlogit}
\alias{logit}
\title{The logit and inverse-logit functions}
\description{
The logit and inverse-logit (also called the logistic function) are
provided.
}
\usage{
invlogit(x)
logit(p)
}
\arguments{
     \item{x}{This is a vector of real values that will be transformed
       to the interval [0,1].}
     \item{p}{This is a vector of probabilities p in the interval [0,1]
       that will be transformed to the real line.}
}
\details{
The \code{logit} function is the inverse of the sigmoid or logistic
function, and transforms a continuous value (usually probability
\eqn{p}) in the interval [0,1] to the real line (where it usually is
the logarithm of the odds). The \code{logit} function is
\eqn{\log(p / (1-p))}{log(p / (1 - p))}.

The \code{invlogit} function (called either the inverse logit or the
logistic function) transforms a real number (usually the logarithm of
the odds) to a value (usually probability \eqn{p}) in the interval
[0,1]. The \code{invlogit} function is
\eqn{\frac{1}{1 + \exp(-x)}}{1 / (1 + exp(-x))}.

If \eqn{p} is a probability, then \eqn{\frac{p}{1-p}}{p/(1-p)} is the
corresponding odds, while the \code{logit} of \eqn{p} is the logarithm
of the odds. The difference between the logits of two probabilities is
the logarithm of the odds ratio. The derivative of probability \eqn{p}
in a logistic function (such as \code{invlogit}) is:
\eqn{\frac{d}{dx} = p(1-p)}{(d / dx) = p * (1 - p)}. 
}
\value{
  \code{invlogit} returns probability \code{p}, and
  \code{logit} returns \code{x}.
}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
x <- -5:5
p <- invlogit(x)
x <- logit(p)
}
\keyword{inverse logit, logistic function, logit}