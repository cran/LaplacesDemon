\name{KLD}
\alias{KLD}
\title{Kullback-Leibler Divergence (KLD)}
\description{
This function calculates the Kullback-Leibler divergence (KLD) between
two probability distributions, and has many uses, especially in prior
elicitation, reference priors, and posterior predictive checks.
}
\usage{
KLD(px, py, base)
}
\arguments{
     \item{px}{This is a required vector of probability densities,
          considered as \eqn{p(\textbf{x})}{p(x)}.}
     \item{py}{This is a required vector of probability densities,
          considered as \eqn{p(\textbf{y})}{p(y)}.}
     \item{base}{This optional argument specifies the logarithmic base,
       which defaults to \code{base=exp(1)} (or \eqn{e}) and represents
       information in nats, where \code{base=2} represents information
       in bits.}
}
\details{
The Kullback-Leibler divergence (KLD) is known by many names, some of
which are Kullback-Leibler distance and K-L. KLD is an asymmetric measure
of the difference, distance, or direct divergence between two
probability distributions \eqn{p(\textbf{y})}{p(y)} and
\eqn{p(\textbf{x})}{p(x)} (Kullback and Leibler, 1951). Here,
\eqn{p(\textbf{y})}{p(y)} represents the ``true'' distribution of data,
observations, or theoretical distribution, and \eqn{p(\textbf{x})}{p(x)}
represents a theory, model, or approximation of \eqn{p(\textbf{y})}{p(y)}.

For probability distributions \eqn{p(\textbf{y})}{p(y)} and
\eqn{p(\textbf{x})}{p(x)} that are discrete (whether the underlying
distribution is continuous or discrete, the observations themselves
are always discrete, such as from \eqn{i=1,\dots,N}{i=1,...,N}),

\deqn{\mathrm{KLD}[p(\textbf{y}) || p(\textbf{x})] = \sum^N_i
  p(\textbf{y}_i) \log\frac{p(\textbf{y}_i)}{p(\textbf{x}_i)}}{KLD[p(y)||p(x)] = sum of p(y[i]) log(p(y[i]) / p(x[i]))}

In Bayesian inference, KLD can be used as a measure of the information
gain in moving from a prior distribution, \eqn{p(\theta)}{p(theta)}, to
a posterior distribution, \eqn{p(\theta | \textbf{y})}{p(theta | y)}. As
such, KLD is the basis of reference priors (Berger, Bernardo, and Sun, 2009).
}
\value{
     \code{KLD} returns a list with the following components:
     \item{KLD.px.py}{This is \eqn{\mathrm{KLD}_i[p(\textbf{x}_i) || p(\textbf{y}_i)]}{KLD[i](p(x[i]) || p(y[i]))}.}
     \item{KLD.py.px}{This is \eqn{\mathrm{KLD}_i[p(\textbf{y}_i) || p(\textbf{x}_i)]}{KLD[i](p(y[i]) || p(x[i]))}.}
     \item{mean.KLD}{This is the mean of the two components above.}
     \item{sum.KLD.px.py}{This is \eqn{\mathrm{KLD}[p(\textbf{x}) || p(\textbf{y})]}{KLD(p(x) || p(y))}.}
     \item{sum.KLD.py.px}{This is \eqn{\mathrm{KLD}[p(\textbf{y}) || p(\textbf{x})]}{KLD(p(y) || p(x))}.}
     \item{mean.sum.KLD}{This is the mean of the two components above.}
     \item{intrinsic.discrepancy}{This is the minimum of the two summed components.}
}
\references{
  Berger, J.O., Bernardo, J.M., and Sun, D. (2009). "The Formal
  Definition of Reference Priors". The Annals of Statistics, 37(2),
  p. 905--938.

  Hall, B. (2012). "Bayesian Inference", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
  
  Kullback, S. and Leibler, R.A. (1951). "On Information and
  Sufficiency". The Annals of Mathematical Statistics, 22(1), p. 79--86.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\examples{
px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
KLD(px,py)
}
\keyword{distribution, elicitation, reference priors}