\name{IVO}
\alias{IVO}
\title{Initial Values Optimization}
\description{
  The Initial Values Optimization \code{IVO} function is a
  deterministic algorithm that maximizes the logarithm of the
  unnormalized joint posterior density with an adaptive gradient ascent
  algorithm. This function is useful for optimizing initial values for
  \code{LaplacesDemon}.
}
\usage{
IVO( Model=NULL, parm=NULL, Data=NULL, Interval=1.0E-6,
     Iterations=100, Stop.Tolerance=1.0E-5)
}
\arguments{
     \item{Model}{
          This required argument receives the model from a user-defined
          function. The user-defined function is where the model is
	  specified. \code{LaplacesDemon} passes two arguments to the
	  model function, \code{parms} and \code{Data}. For more
	  information, see the \code{LaplacesDemon} function and
	  vignette.
          }
     \item{parm}{
          This argument requires a vector of initial values equal in
          length to the number of parameters. \code{IVO} will attempt to
	  optimize these initial values for the parameters, for later
	  use with the \code{LaplacesDemon} function.
          }
     \item{Data}{
          This required argument accepts a list of data.
          }
     \item{Interval}{
          This argument receives an interval for estimating approximate
          gradients. The logarithm of the unnormalized joint posterior
	  density of the Bayesian model is evaluated at +/- 1/2 of this
	  interval. \code{Interval} defaults to 1.0E-6.
          }
     \item{Iterations}{
          This argument accepts an integer that determines the number
          of iterations that \code{IVO} will attempt to maximize the
	  logarithm of the unnormalized joint posterior density.
	  \code{Iterations} defaults to 100. \code{IVO} will stop before
	  this number of iterations if the tolerance is less than or equal
	  to the \code{Stop.Tolerance} criterion (see below).
          }
     \item{Stop.Tolerance}{
	  This argument accepts any positive number, and defaults to 
          1.0E-5. At each iteration, the square root of the sum of the
	  squared differences of the logarithm of the unnormalized joint
	  posterior density is calculated. If this result is less than
	  or equal to the value of \code{Stop.Tolerance}, then
	  \code{IVO} has converged to the user-specified tolerance and
	  will terminate at the end of the current iteration. 
	  }
}
\details{
     \code{IVO} may be called by the user before using
     \code{LaplacesDemon}, or \code{LaplacesDemon} may call this
     function if all initial values are zero. Chasing convergence with
     \code{IVO} may be time-consuming and unimportant. The goal,
     instead, is to improve the logarithm of the unnormalized joint
     posterior density so that it is easier for the \code{LaplacesDemon}
     function to begin updating the parameters in search of the target
     distributions. This can be difficult when the initial values are in
     low-probability regions, and can cause unreasonably low acceptance
     rates.

     \code{IVO} uses a first-order optimization algorithm that seeks a
     global maximum of the logarithm of the unnormalized joint posterior
     density by taking adaptive steps proportional to an approximate
     gradient. This is a gradient ascent algorithm, and is called a
     gradient descent or steepest descent algorithm elsewhere for
     minimization. Laplace's Demon uses the \code{IVO} algorithm to
     optimize initial values and save time for the user.

     This algorithm assumes that the logarithm of the unnormalized joint
     posterior density is defined and differentiable. An approximate
     gradient is taken for each initial value as the difference in the
     logarithm of the unnormalized joint posterior density due to a
     slight increase versus decrease in the parameter.

     At 10 evenly-spaced times, \code{IVO} attempts several step sizes,
     which are also called rate parameters in other literature, and
     selects the best step size. Thereafter, each iteration in which an
     improvement does not occur, the step size shrinks, being multiplied
     by 0.999.

     An nth order approximation has an order of magnitude of error of
     x^n, where a Taylor polynomial of n degrees is an nth order
     approximation. SAGA is only a first-order, or linear,
     approximation, and is equivalent to a polynomial of degree one. For
     example, Newton-Rhapson is a second-order algorithm, using
     quadratic polynomials due to the Hessian matrix.

     Gradient ascent is criticized for sometimes being relatively slow
     when close to the maximum, and its asymptotic rate of convergence
     is inferior to other methods. However, compared to other popular
     optimization algorithms such as Newton-Rhapson, an advantage of the
     gradient ascent is that it works in infinite dimensions, requiring
     only sufficient computer memory. Although Newton-Rhapson converges
     in fewer iterations and is a second-order approximation,
     calculating the inverse of the Hessian matrix of second-derivatives
     is more computationally expensive and subject to singularities.
     Therefore, gradient ascent takes longer to converge, but is more
     generalizable.
}
\value{
     \code{IVO} returns a list with the following components:
     \item{Converged}{
          This is a logical indicator of whether or not \code{IVO} 
          converged within the specified \code{Iterations} according to
          the supplied \code{Stop.Tolerance} criterion. Convergence does
          not indicate that the global maximum has been found, but only
	  that the tolerance was less than or equal to the
	  \code{Stop.Tolerance} criterion. 
          }
     \item{LP.Final}{
          This reports the final scalar value for the logarithm of the
          unnormalized joint posterior density.
          }
     \item{LP.Initial}{
          This reports the initial scalar value for the logarithm of the
          unnormalized joint posterior density. 
          }
     \item{Parms.Final}{
          This is the vector of parameters that is ready to passed to
          the \code{LaplacesDemon} function as initial values.
	  }
     \item{Parms.Initial}{
	  This is the vector of initial values that were originally 
          given to \code{IVO} in the \code{parm} argument.
	  }
     \item{Step.Size.Final}{
          This is the final, scalar \code{Step.Size} value at the end 
          of the \code{IVO} algorithm.
	  }
     \item{Step.Size.Initial}{
          This is the initial, scalar \code{Step.Size}.
          }
     \item{Stop.Tolerance}{
          This is the \code{Stop.Tolerance} criterion.
	  }
     \item{Tolerance}{
	  This is the last \code{Tolerance} of the \code{IVO}
	  algorithm. It is calculated as the square root of the sum of
	  the squared differences between a new and current vector of
	  parameters.
          }
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
N <- NROW(demonsnacks)
J <- NCOL(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(1,3:10)]))
for (j in 2:J) {X[,j] <- (X[,j] - mean(X[,j])) / (2*sd(X[,j]))}
parm.names <- rep(NA,J+1)
for (j in 1:J) {parm.names[j] <- paste("beta[",j,"]",sep="")}
parm.names[J+1] <- "log.tau"
MyData <- list(J=J, X=X, parm.names=parm.names, y=y)

############################  Initial Values  #############################
Initial.Values <- c(rep(0,J), log(1))

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Prior Parameters
     beta.mu <- rep(0,Data$J)
     beta.tau <- rep(1.0E-3,Data$J)
     tau.alpha <- 1.0E-3
     tau.beta <- 1.0E-3
     ### Parameters
     beta <- parm[1:Data$J]
     tau <- exp(parm[Data$J+1])
     ### Log of Prior Densities
     beta.prior <- dnorm(beta, beta.mu, 1/sqrt(beta.tau), log=TRUE)
     tau.prior <- dgamma(tau, tau.alpha, tau.beta, log=TRUE)
     ### Log-Posterior
     mu <- beta \%*\% t(Data$X)
     LL <- sum(dnorm(Data$y, mu, 1/sqrt(tau), log=TRUE))
     LP <- LL + sum(beta.prior) + tau.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(tau,mu[1]), yhat=mu,
          parm=parm)
     return(Modelout)
     }

Fit <- IVO(Model, Initial.Values, Data=MyData)
Fit
Fit <- IVO(Model, Initial.Values, Data=MyData, Iterations=200)
Fit
#Fit$Parms.Final are submitted to LaplacesDemon as Initial.Values.
#End
}
\keyword{
Adaptive, Bayesian Inference, Gradient Ascent, Initial Values, IVO,
LaplacesDemon, Laplace's Demon, Optimization
}
