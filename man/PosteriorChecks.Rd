\name{PosteriorChecks}
\alias{PosteriorChecks}
\title{Posterior Checks}
\description{
Not to be confused with posterior predictive checks, this function
provides some posterior checks, such as the probability that each theta
is greater than zero, the kurtosis and skewness of the posterior
distributions, and a correlation matrix is provided.
}
\usage{
PosteriorChecks(x, Parms)
}
\arguments{
    \item{x}{This required argument accepts an object of class demonoid
      or laplace.}
    \item{Parms}{This argument accepts a vector of quoted strings to be
      matched for selecting parameters. This argument defaults to
      \code{NULL} and selects every parameter. Each quoted string is
      matched to one or more parameter names with the \code{grep}
      function. For example, if the user specifies
      \code{Parms=c("eta","tau")}, and if the parameter names are
      beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
      will be selected, because the string \code{eta} is within
      \code{beta}. Since \code{grep} is used, string matching uses
      regular expressions, so beware of meta-characters, though these
      are acceptable: ".", "[", and "]".}
}
\details{\code{PosteriorChecks} is a supplemental function that returns
  a list with two components.}
\value{
     \code{Posterior.Correlation} returns a correlation matrix of the
     parameters selected with the \code{Parms} argument.
     \code{Posterior.Summary} returns a matrix in which each row is a
     parameter and there are three columns: p(theta > 0), Kurtosis, and
     Skewness. The first column, p(theta > 0), indicates the importance
     of a parameter by reporting how much of the distribution overlaps
     zero. Kurtosis and skewness are useful posterior checks that may
     suggest that a posterior distribution is not normal, for example.
}
\keyword{posterior, check}