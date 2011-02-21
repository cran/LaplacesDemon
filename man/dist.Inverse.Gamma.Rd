\name{dist.Inverse.Gamma}
\alias{dinvgamma}
\alias{rinvgamma}
\title{The Inverse Gamma Distribution}
\description{
  This is the density function and random generation from the inverse
  gamma distribution.
}
\usage{
dinvgamma(x, shape, scale=1, log=FALSE)
rinvgamma(n, shape, scale=1)
}
\arguments{
  \item{n}{This is the number of draws from the distribution.}
  \item{x}{This is the scalar location to evaluate density.}
  \item{shape}{This is the scalar shape parameter.}
  \item{scale}{This is the scalar scale parameter (which defaults to
    one).}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\value{
  \code{dinvgamma} evaluates the density at \code{x}. \code{rinvgamma} takes
  \code{n} draws from the inverse Gamma distribution. The parameterization is
  consistent with the Gamma Distribution in the stats package.
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: p(theta) = (beta^alpha / gamma(alpha)) * theta^(-(alpha + 1)) *exp(-beta / theta), theta > 0
    \item Inventor: Unknown (to me, anyway)
    \item Notation 1: theta ~ Inv-gamma(alpha, beta)
    \item Notation 2: p(theta) = Inv-gamma(theta | alpha, beta)
    \item Parameter 1: shape alpha > 0
    \item Parameter 2: scale beta > 0
    \item Mean: E(theta) = beta / (alpha - 1), for alpha > 1
    \item Variance: var(theta) = beta^2 / ((alpha - 1)^2 * (alpha - 2)), alpha > 2
    \item Mode: mode(theta) = beta / (alpha + 1)
}

The inverse-gamma is the conjugate prior distribution for the normal
or Gaussian variance. The density is always finite; its integral is
finite if \code{alpha > 0}. An uninformative distribution is obtained as
\code{alpha, beta -> 0}.

These functions are similar to those in the \code{MCMCpack} package.
}
\examples{
x <- dinvgamma(4.3, 1.1)
x <- rinvgamma(10, 3.3)
}
\keyword{distribution}

