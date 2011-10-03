\name{PosteriorChecks}
\alias{PosteriorChecks}
\title{Posterior Checks}
\description{
Not to be confused with posterior predictive checks, this function
provides some posterior checks, such as the probability that each
posterior parameter (referred to generically as \eqn{\theta}{theta}),
is greater than zero [\eqn{p(\theta > 0)}{p(theta > 0)}], the kurtosis
and skewness of the posterior distributions, and a correlation matrix
is provided.
}
\usage{
PosteriorChecks(x, Parms)
}
\arguments{
    \item{x}{This required argument accepts an object of class
      \code{demonoid} or \code{laplace}.}
    \item{Parms}{This argument accepts a vector of quoted strings to be
      matched for selecting parameters. This argument defaults to
      \code{NULL} and selects every parameter. Each quoted string is
      matched to one or more parameter names with the \code{grep}
      function. For example, if the user specifies
      \code{Parms=c("eta", "tau")}, and if the parameter names are
      beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
      will be selected, because the string \code{eta} is within
      \code{beta}. Since \code{grep} is used, string matching uses
      regular expressions, so beware of meta-characters, though these
      are acceptable: ".", "[", and "]".}
}
\details{\code{PosteriorChecks} is a supplemental function that returns
  a list with two components.}
\value{
     \code{PosteriorChecks} returns an object that is a list with the
          following components:
     \item{Posterior.Correlation}{
	  This is a correlation matrix of the parameters selected with
	  the \code{Parm} argument.
	  }
     \item{Posterior.Summary}{This is a matrix in which each row is a
          parameter and there are three columns: p(theta > 0), Kurtosis,
	  and Skewness. The first column, p(theta > 0), indicates the
	  importance of a parameter by reporting how much of the
	  distribution is greater than zero. An important distribution
	  will have a result at least as extreme as 0.025 or 0.975, and
	  an unimportant distribution is centered at 0.5. Kurtosis and
	  skewness are useful posterior checks that may suggest that a
	  posterior distribution is not normal or does not fit well with
	  a distributional assumption, assuming a distributional
	  assumption exists, which it may not. 
          }
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{### See the LaplacesDemon function for an example.}
\keyword{posterior, check}