\name{LaplaceApproximation}
\alias{LaplaceApproximation}
\title{Laplace Approximation}
\description{
  The \code{LaplaceApproximation} function deterministically maximizes
  the logarithm of the unnormalized joint posterior density with an
  adaptive gradient ascent algorithm. The goal of Laplace Approximation
  is to estimate the posterior mode and variance of each parameter.
  This function is useful for optimizing initial values and estimating
  a covariance matrix to be input into the \code{LaplacesDemon}
  function, or sometimes for model estimation in its own right.
}
\usage{
LaplaceApproximation(Model=NULL, parm=NULL, Data=NULL, Interval=1.0E-6,
     Iterations=100, Stop.Tolerance=1.0E-5)
}
\arguments{
     \item{Model}{
          This required argument receives the model from a user-defined
          function. The user-defined function is where the model is
	  specified. \code{LaplaceApproximation} passes two arguments to
	  the model function, \code{parms} and \code{Data}. For more
	  information, see the \code{LaplacesDemon} function and
	  vignette.
          }
     \item{parm}{
          This argument requires a vector of initial values equal in
          length to the number of parameters. \code{LaplaceApproximation} 
	  will attempt to optimize these initial values for the
	  parameters, where the optimized values are the posterior
	  modes, for later use with the \code{LaplacesDemon} function.
          }
     \item{Data}{
          This required argument accepts a list of data. The list of
	  data must include \code{mon.names} which contains monitored
	  variable names, and \code{parm.names} which contains parameter
	  names. \code{LaplaceApproximation} must be able to determine
	  the sample size of the data, and will look for a scalar sample
	  size variable \code{n} or \code{N}. If not found, it will look
	  for variable \code{y} or \code{Y}, and attempt to take its
	  number of rows as sample size. \code{LaplaceApproximation}
	  needs to determine sample size due to the asymptotic nature
	  of this method.
          }
     \item{Interval}{
          This argument receives an interval for estimating approximate
          gradients. The logarithm of the unnormalized joint posterior
	  density of the Bayesian model is evaluated at +/- 1/2 of this
	  interval. \code{Interval} defaults to 1.0E-6.
          }
     \item{Iterations}{
          This argument accepts an integer that determines the number
          of iterations that \code{LaplaceApproximation} will attempt to
	  maximize the logarithm of the unnormalized joint posterior density.
	  \code{Iterations} defaults to 100. \code{LaplaceApproximation}
	  will stop before this number of iterations if the tolerance is
	  less than or equal to the \code{Stop.Tolerance} criterion
	  (see below).
          }
     \item{Stop.Tolerance}{
	  This argument accepts any positive number, and defaults to 
          1.0E-5. At each iteration, the square root of the sum of the
	  squared differences of the logarithm of the unnormalized joint
	  posterior density is calculated. If this result is less than
	  or equal to the value of \code{Stop.Tolerance}, then
	  \code{LaplaceApproximation} has converged to the user-specified 
	  tolerance and will terminate at the end of the current iteration. 
	  }
}
\details{
  The Laplace Approximation or Laplace Method is a family of asymptotic
  techniques used to approximate integrals. Laplace's method seems to
  accurately approximate uni-modal posterior moments and marginal
  posterior distributions in many cases. Since it is not applicable in
  all cases, it is recommended here that Laplace Approximation is used
  cautiously in its own right, or preferably, it is used before MCMC.

  After introducing the Laplace Approximation (Laplace, 1774,
  p. 366--367), a proof was published later (Laplace, 1814) as part of
  a mathematical system of inductive reasoning based on probability.
  Laplace used this method to approximate posterior moments.

  Since its introduction, the Laplace Approximation has been applied
  successfully in many disciplines. In the 1980s, the Laplace
  Approximation experienced renewed interest, especially in statistics,
  and some improvements in its implementation were introduced (Tierney
  et al., 1986; Tierney et al., 1989). Only since the 1980s has the
  Laplace Approximation been seriously considered by statisticians in
  practical applications.

  There are many variations of Laplace Approximation, with an effort
  toward replacing Markov chain Monte Carlo (MCMC) algorithms as the
  dominant form of numerical approximation in Bayesian inference. The
  run-time of Laplace Approximation is a little longer than Maximum
  Likelihood Estimation (MLE), and much shorter than MCMC (Azevedo and
  Shachter, 1994). In the \code{LaplacesDemon} package, Laplace
  Approximation may iterate faster or slower than MCMC, so this is not
  the fastest possible implementation of Laplace Approximation. Laplace
  Approximation extends MLE, but shares similar limitations, such as
  its asymptotic nature with respect to sample size. Bernardo and Smith
  (2000) note that Laplace Approximation is an attractive numerical
  approximation algorithm, and will continue to develop, though it
  currently works best with few parameters.

  The \code{LaplaceApproximation} function may be called by the user
  before using \code{LaplacesDemon}, or \code{LaplacesDemon} may call
  this function if all initial values are zero and the sample size is at
  least five times the number of parameters. Chasing convergence with
  \code{LaplaceApproximation} may be time-consuming and unimportant. The
  goal, instead, is to improve the logarithm of the unnormalized joint
  posterior density so that it is easier for the \code{LaplacesDemon}
  function to begin updating the parameters in search of the target
  distributions. This can be difficult when the initial values are in
  low-probability regions, and can cause unreasonably low acceptance rates.

  \code{LaplaceApproximation} seeks a global maximum of the logarithm of
  the unnormalized joint posterior density by taking adaptive steps
  proportional to an approximate gradient. This portion of the
  \code{LaplaceApproximation} function uses a gradient ascent algorithm,
  and is called a gradient descent or steepest descent algorithm
  elsewhere for minimization problems. Laplace's Demon uses the
  \code{LaplaceApproximation} algorithm to optimize initial values,
  estimate covariance, and save time for the user.

  This algorithm assumes that the logarithm of the unnormalized joint
  posterior density is defined and differentiable. An approximate
  gradient is taken for each initial value as the difference in the
  logarithm of the unnormalized joint posterior density due to a slight
  increase versus decrease in the parameter.

  At 10 evenly-space times, \code{LaplaceApproximation} attempts several
  step sizes, which are also called rate parameters in other literature,
  and selects the best step size from a set of 10 fixed options.
  Thereafter, each iteration in which an improvement does not occur, the
  step size shrinks, being multiplied by 0.999.

  Gradient ascent is criticized for sometimes being relatively slow when
  close to the maximum, and its asymptotic rate of convergence is
  inferior to other methods. However, compared to other popular
  optimization algorithms such as Newton-Rhapson, an advantage of the
  gradient ascent is that it works in infinite dimensions, requiring
  only sufficient computer memory. Although Newton-Rhapson converges in
  fewer iterations, calculating the inverse of the negative Hessian
  matrix of second-derivatives is more computationally expensive and
  subject to singularities. Therefore, gradient ascent takes longer to
  converge, but is more generalizable.

  After \code{LaplaceApproximation} finishes, due either to early
  convergence or completing the number of specified iterations, it
  approximates the Hessian matrix of second derivatives, and attempts to
  calculate the covariance matrix by taking the inverse of the negative
  of this matrix. If successful, then this covariance matrix may be
  passed to \code{LaplacesDemon}, and the diagonal of this matrix is the
  variance of the parameters. If unsuccessful, a scaled identity matrix
  is returned, and each parameter's variance will be 1.

}
\value{
     \code{LaplaceApproximation} returns an object of class \code{laplace}
          that is a list with the following components:
     \item{Call}{
	  This is the matched call of \code{LaplaceApproximation}.
	  }
     \item{Converged}{
          This is a logical indicator of whether or not
          \code{LaplaceApproximation} converged within the specified
	  \code{Iterations} according to the supplied
	  \code{Stop.Tolerance} criterion. Convergence does not indicate
	  that the global maximum has been found, but only that the
	  tolerance was less than or equal to the \code{Stop.Tolerance}
	  criterion. 
          }
     \item{Covar}{
          This covariance matrix is the negative inverse of an
	  approximate Hessian matrix. If an error is encountered in
	  attempting to solve the matrix inversion, then an
	  identity matrix is returned. The \code{Covar} matrix may be
	  scaled and input into the \code{Covar} argument of the
	  \code{LaplacesDemon} function for further MCMC estimation, or
	  the diagonal of this matrix may be used to represent the
	  posterior variance of the parameters, provided the algorithm
	  converged and matrix inversion was successful. To scale this
	  matrix for use with Laplace's Demon, multiply it by
	  \code{2.38^2/d}, where \code{d} is the number of initial
	  values.
          }
     \item{Deviance}{
	  This is a vector of the iterative history of the deviance in
          the \code{LaplaceApproximation} function, as it sought
	  convergence. 
	  }
     \item{History}{
	  This is a matrix of the iterative history of the parameters in
          the \code{LaplaceApproximation} function, as it sought
	  convergence. 
	  }
     \item{Initial.Values}{
	  This is the vector of initial values that was originally 
          given to \code{LaplaceApproximation} in the \code{parm}
	  argument.
	  }
     \item{LML}{
          This is an approximation of the logarithm of the marginal
          likelihood of the data, and is the logarithmic form of
	  equation 4 in Lewis and Raftery (1997). As a rough estimate of
	  Kass and Raftery (1995), \code{LML} is worrisome when the
	  sample size of the data is less than five times the number of
	  parameters, and \code{LML} should be adequate in most problems
	  when the sample size of the data exceeds twenty times the
	  number of parameters (p. 778). \code{LML} is useful for
	  comparing multiple models with the \code{BayesFactor}
	  function. The Laplace-Metropolis estimator is inappropriate
	  with hierarchical models, and the extension (the Compound
	  Laplace-Metropolis estimator) will soon be attempted to be
	  included.
          }
     \item{LP.Final}{
          This reports the final scalar value for the logarithm of the
          unnormalized joint posterior density.
          }
     \item{LP.Initial}{
          This reports the initial scalar value for the logarithm of the
          unnormalized joint posterior density. 
          }
     \item{Minutes}{
          This indicates the number of minutes that
          \code{LaplaceApproximation} was running, and this includes the
	  initial checks as well as creating summaries.
	  }
     \item{Step.Size.Final}{
          This is the final, scalar \code{Step.Size} value at the end 
          of the \code{LaplaceApproximation} algorithm.
	  }
     \item{Step.Size.Initial}{
          This is the initial, scalar \code{Step.Size}.
          }
     \item{Stop.Tolerance}{
          This is the \code{Stop.Tolerance} criterion.
	  }
     \item{Summary}{
          This is a summary matrix. Rows are parameters. The following
          columns are included: Mode, SD (Standard Deviation), LB
	  (Lower Bound), and UB (Upper Bound).
	  }
     \item{Tolerance}{
	  This is the last \code{Tolerance} of the \code{LaplaceApproxiation}
	  algorithm. It is calculated as the square root of the sum of
	  the squared differences between a new and current vector of
	  parameters.
          }
}
\references{
     Azevedo-Filho, A. and Shachter, R. (1994). "Laplace's Method
     Approximations for Probabilistic Inference in Belief Networks with
     Continuous Variables". In "Uncertainty in Artificial Intelligence",
     Mantaras, R. and Poole, D., Morgan Kauffman, San Francisco, CA,
     p. 28--36.
     
     Bernardo, J.M. and Smith, A.F.M. (2000). "Bayesian Theory". John
     Wiley \& Sons: West Sussex, England.

     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}

     Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". Journal of the
     American Statistical Association, 90(430), p. 773--795.
     
     Laplace, P. (1774). "Memoire sur la Probabilite des Causes par les
     Evenements." l'Academie Royale des Sciences, 6, 621--656. English
     translation by S.M. Stigler in 1986 as "Memoir on the Probability
     of the Causes of Events" in Statistical Science, 1(3), 359--378.

     Laplace, P. (1814). "Essai Philosophique sur les Probabilites."
     English translation in Truscott, F.W. and Emory, F.L. (2007) from
     (1902) as "A Philosophical Essay on Probabilities". ISBN
     1602063281, translated from the French 6th ed. (1840).

     Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
     Posterior Simulation with the Laplace-Metropolis
     Estimator". Journal of the American Statistical Association, 92,
     p. 648--655.
     
     Tierney, L. and Kadane, J.B. (1986). "Accurate Approximations for
     Posterior Moments and Marginal Densities". Journal of the American
     Statistical Association, 81(393), p. 82--86.

     Tierney, L., Kass. R., and Kadane, J.B. (1989). "Fully Exponential
     Laplace Approximations to Expectations and Variances of Nonpositive
     Functions". Journal of the American Statistical Association,
     84(407), p. 710--716.
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
J <- 6
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(1,3:6)]))
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("tau","mu[1]")
parm.names <- parm.names(list(beta=rep(0,J), log.tau=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)

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

Fit <- LaplaceApproximation(Model, Initial.Values, Data=MyData, Iterations=1000)
Fit
print(Fit)
caterpillar.plot(Fit, Parms=1:6)
plot(Fit,MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData)
summary(Pred)
plot(Pred)
#Fit$Covar is scaled (2.38^2/d) and submitted to LaplacesDemon as Covar.
#Fit$Summary[,1] is submitted to LaplacesDemon as Initial.Values.
#End
}
\keyword{
Adaptive, Bayesian Inference, Gradient Ascent, Initial Values,
Laplace Approximation, LaplaceApproximation, LaplacesDemon,
Laplace's Demon, Optimization
}
