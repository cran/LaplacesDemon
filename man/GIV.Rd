\name{GIV}
\alias{GIV}
\title{Generate Initial Values}
\description{
  The \code{GIV} function generates initial values for use with the
  \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}}
  functions.
}
\usage{
GIV(Model, Data, n=1000)
}
\arguments{
  \item{Model}{This required argument is a model specification
    function. For more information, see \code{\link{LaplacesDemon}}.}
  \item{Data}{This required argument is a list of data. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{n}{This is the number of attempts to generate acceptable
    initial values.}
}
\details{
  Initial values are required for optimization algorithms. A user may
  specify initial values, or use the \code{GIV} function for random
  generation. Initial values determined by the user may fail to produce
  a finite posterior in complicated models, and the \code{GIV} function
  is here to help.

  \code{GIV} has several uses. First, the
  \code{\link{LaplaceApproximation}} and \code{\link{LaplacesDemon}}
  functions use \code{GIV} internally if unacceptable initial values are
  discovered. Second, the user may use \code{GIV} when developing their
  model specification function, \code{Model}, to check for potential
  problems. Third, the user may prefer to randomly generate acceptable
  initial values. Lastly, \code{GIV} is recommended when running
  multiple chains (such as for later use with the
  \code{Gelman.Diagnostic}) for dispersed starting locations.

  In \code{GIV}, initial values are constrained to the interval
  [-100,100]. This is done to prevent numeric overflows with parameters
  that are exponentiated within the model specification function. First,
  \code{GIV} passes the upper and lower bounds of this interval to the
  model, and any changed parameters are noted. Initial values are
  randomly generated with a normal proposal and a small variance at the
  center of the returned range of each parameter. As \code{GIV} fails to
  find acceptable initial values, the algorithm iterates toward its
  maximum number of iterations, \code{n}. In each iteration, the
  variance increases for the proposal.
  
  Initial values are considered acceptable only when the first two
  returned components of \code{Model} (which are \code{LP} and
  \code{Dev}) are finite.

  If \code{GIV} fails to return acceptable initial values, then it is
  best to study the model specification function. When the model is
  complicated, here is a suggestion. Remove the log-likelihood,
  \code{LL}, from the equation that calculates the logarithm of the
  unnormalized joint posterior density, \code{LP}. For example, convert
  \code{LP <- LL + beta.prior} to \code{LP <- beta.prior}. Now, maximize
  \code{LP}, which is merely the set of prior densities, with any
  optimization algorithm. Replace \code{LL}, and run the model with
  initial values that are in regions of high prior density. If this
  fails, then the model specification should be studied closely, because
  a non-finite posterior should (especially) never be associated with
  regions of high prior density.
}
\value{
  The \code{GIV} function returns a vector equal in length to the
  number of parameters, and each element is an initial value for the
  associated parameter in \code{Data$parm.names}. When \code{GIV} fails
  to find acceptable initial values, each returned element is \code{NA}.
}
\references{
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{as.initial.values}},
  \code{\link{Gelman.Diagnostic}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
N <- nrow(demonsnacks)
J <- ncol(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(1,3:10)]))
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[1:Data$J]
     sigma <- exp(parm[Data$J+1])
     ### Log of Prior Densities
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Likelihood
     mu <- tcrossprod(beta, Data$X)
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma), yhat=mu,
          parm=parm)
     return(Modelout)
     }

########################  Generate Initial Values  ########################
Initial.Values <- GIV(Model, MyData)
}
\keyword{utility}
