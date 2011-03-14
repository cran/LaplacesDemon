\name{LaplacesDemon}
\alias{LaplacesDemon}
\title{Laplace's Demon}
\description{
The \code{LaplacesDemon} function is the main function of Laplace's
Demon. Given data, a model specification, and initial values,
\code{LaplacesDemon} maximizes the logarithm of the unnormalized joint
posterior density and provides samples of the marginal posterior
distributions, deviance, and other monitored variables.
}
\usage{
LaplacesDemon(Model=NULL, Data=NULL, Adaptive=0, Covar=NULL, DR=1,
     Initial.Values=NULL, Iterations=1000, Periodicity=1, Status=100,
     Thinning=1)
}
\arguments{
     \item{Model}{
          This required argument receives the model from a user-defined
          function. The user-defined function is where the model is
	  specified. \code{LaplacesDemon} passes two arguments to the
	  model function, \code{parms} and \code{Data}, and receives
	  five arguments from the model function: \code{LP} (the logarithm of
	  the unnormalized joint posterior), \code{Dev} (the deviance),
	  \code{Monitor} (the monitored variables), \code{yhat} (the
	  variables for posterior predictive checks), and \code{parm},
	  the vector of parameters, which may be constrained in the
	  model function.
          }
     \item{Data}{
          This required argument accepts a list of data. The list of
	  data must contain \code{mon.names} which contains monitored
	  variable names, and must contain \code{parm.names} which
	  contains parameter names.
          } 
     \item{Adaptive}{
          This argument instructs Laplace's Demon to adapt the proposal
          covariance matrix, or proposal variance, based on the history of
          each chain, beginning with a specified iteration. If \code{Adaptive=0}
          or is greater than the number of iterations, then adaptation will 
          not occur. For example, if \code{Adaptive=100}, then the proposal
          variances may begin to adapt at the 100th iteration. If the
          argument is unspecified, then it defaults to \code{Adaptive=0}, and does 
          not adapt. Note that adaptation also depends on the \code{Periodicity}
	  argument (see below).
          }
     \item{Covar}{
          This argument defaults to \code{NULL}, but may otherwise accept a K x K
          proposal covariance matrix for the first adaptation of the proposal 
          covariances, where K is the number of dimensions (or parameters). 
          The proposal covariance matrix will be re-estimated with each
          adaptation according to the entire history of all chains. When
          the model is updated for the first time, \code{NULL} should be used,
          unless there is a better estimate at the variance of each
          target distribution as well as the associated covariances. Once
          Laplace's Demon has finished updating, it may be desired to
          continue updating where it left off, in which case the proposal
          covariance matrix from the last run can be input into the next
	  run. The covariance matrix may also be input from the
	  \code{LaplaceApproximation} function, if used.
        }
     \item{DR}{
	  This argument indicates whether or not Delayed Rejection will
          occur, and it defaults to \code{DR=1}. If Delayed Rejection
          does occur, then a rejected proposal will be followed by
          another proposal made with half the covariance. Otherwise,
	  rejections are handled in the traditional way for a Metropolis
	  sampler.
	  }
     \item{Initial.Values}{
          This argument requires a vector of initial values equal in
          length to the number of parameters. Each initial value will be
          the starting point for an adaptive chain or a non-adaptive
	  Markov chain of a parameter. If all initial values are set to
	  zero, then Laplace's Demon will attempt to optimize the
	  initial values with the (\code{LaplaceApproximation}) 
	  function. After Laplace's Demon finishes updating, it may be
	  desired to continue updating from where it left off. To
	  continue, this argument should receive the last iteration of
	  the previous update. For example, if Fit is the output
	  object, then 
          \code{Initial.Values=Fit$Posterior1[NROW(Fit$Posterior1),]}.
          }
     \item{Iterations}{
          This required argument accepts integers greater than 10 and up 
          to 1,000,000, and determines the number of iterations that
	  Laplace's Demon will update the parameters while searching for
	  target distributions.
          }
     \item{Periodicity}{
          This argument accepts positive integers and determines the
          periodicity of adaptation in terms of iterations. For example,
          if \code{Periodicity=10}, then adaptation will occur every 10
          iterations, and begins after the number of iterations specified
          in the \code{Adaptive} argument (above). When
	  \code{Periodicity=1}, it will adapt every iteration. Frequent
	  adaptation is desirable, but requires more time due to
	  calculations.
          }
     \item{Status}{
          This argument accepts integers between 1 and the number of 
          iterations, and indicates how often the user would like the status
          of the number of iterations and proposal type (multivariate or 
          single-component) printed to the screen. For example, if a model
          is updated for 1,000 iterations and \code{Status=200}, then a
	  status message will be printed at the following iterations:
	  200, 400, 600, and 800.
          }
     \item{Thinning}{
          This argument accepts integers between 1 and the number of
          iterations, and indicates that every nth iteration will be
          retained, while the other iterations are discarded. If
          \code{Thinning=5}, then every 5th iteration will be retained.
	  Thinning is performed to reduce autocorrelation and the number of
          marginal posterior samples.
	  }
}
\details{
     \code{LaplacesDemon} offers four customized MCMC algorithms for 
     numerical approximation in Bayesian inference. The four algorithms
     are 1) Adaptive Metropolis (AM), 2) Delayed Rejection Metropolis
     (DRM), 3) Delayed Rejection Adaptive Metropolis (DRAM), and 4)
     Random-Walk Metropolis (RWM). References are provided below.
     All algorithms attempt multivariate proposals, but are modified
     to use single-component proposals if multivariate estimation is
     problematic. Additionally, the \code{LaplaceApproximation} function
     may be used to attempt to optimize initial values with a gradient
     ascent algorithm, which the \code{LaplacesDemon} function attempts
     automatically when all initial values are set to zero and sample
     size is at least five times the number of parameters. For
     more information, see the accompanying vignettes entitled 
     "Bayesian Inference" and "LaplacesDemon Tutorial".
}
\value{
     \code{LaplacesDemon} returns an object of class \code{demonoid}
     that is a list with the following components:
     \item{Acceptance.Rate}{
          This is the acceptance rate of the MCMC algorithm, indicating
          percentage of iterations in which the proposals were
          accepted. The optimal acceptance rate varies with the 
          number of parameters, ranging from 0.44 for one parameter (one
	  IID Gaussian target distribution) to 0.234 for an infinite
	  number of parameters (IID Gaussian target distributions), and
	  0.234 is approached quickly as the number of parameters
	  increases.
          }
     \item{Adaptive}{
          This reports the value of the \code{Adaptive} argument.
          }
     \item{Algorithm}{
          This reports the specific algorithm used. 
          }
     \item{Call}{
	  This is the matched call of \code{LaplacesDemon}.
	  }
     \item{Covar}{
	  This stores the d x d proposal covariance matrix of the most recent
          adaptation, where d is the dimension or number of parameters
	  or initial values. If the model is updated in the future, then this 
          matrix can be used to start the next update where the last update
	  left off. Only the diagonal of this matrix is reported in the
	  associated \code{print} function.
	  }
     \item{CovarDHis}{
	  This n x d matrix stores the diagonal of the proposal
	  covariance matrix of each adaptation in each of n rows for d
	  dimensions, where the dimension is the number of parameters or
	  length of the initial values vector. The proposal covariance
	  matrix should change less over time.
	  }
     \item{Deviance}{
          This is a vector of the deviance of the model, with a length
          equal to the number of thinned samples that were retained.
          Deviance is useful for considering model fit, and is equal to
	  the sum of the log-likelihood for all rows in the data set,
	  which is then multiplied by negative two.
	  }
     \item{DIC1}{
	  This is a vector of three values: Dbar, pD, and DIC. Dbar is
          the mean deviance, pD is a measure of model complexity
	  indicating the effective number of parameters, and DIC is the
	  Deviance Information Criterion, which is a model fit statistic
	  that is the sum of Dbar and pD. \code{DIC1} is calculated
	  over all retained samples. Note that pD is calculated as
	  \code{var(Deviance)/2} as in Gelman et al. (2004).
          }
     \item{DIC2}{
          This is identical to \code{DIC1} above, except that it is
          calculated over only the samples that were considered by the
          \code{Geweke.Diagnostic} to be stationary for all parameters. 
	  If stationarity (or a lack of trend) is not estimated for all
	  parameters, then \code{DIC2} is set to missing values.
	  }
     \item{DR}{
	  This reports the value of the \code{DR} argument.
          }
     \item{Initial.Values}{
          This is the vector of \code{Initial.Values}, which may have 
          been optimized with the \code{LaplaceApproximation} function.
          }
     \item{Iterations}{
          This reports the number of \code{Iterations} for updating.
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
     \item{Minutes}{
          This indicates the number of minutes that \code{LaplacesDemon}
          was running, and this includes the initial checks as well as
	  time it took the \code{LaplaceApproximation} function,
	  assessing stationarity, effective sample size (ESS), and
	  creating summaries.
	}
     \item{Model}{
          This contains the model specification \code{Model}.
	  }
     \item{Monitor}{
          This is a vector or matrix of one or more monitored variables,
          which are variables that were specified in the
	  \code{Model} function to be observed as chains (or
	  Markov chains, if \code{Adaptive=0}), but that were not
          deviance or parameters.
	  }
     \item{Parameters}{
          This reports the number of parameters.
          }
     \item{Periodicity}{
          This reports the value of the \code{Periodicity} argument.
          }
     \item{Posterior1}{
          This is a matrix of marginal posterior distributions composed of 
          thinned samples, with a number of rows equal to the number of
	  thinned samples and a number of columns equal to the number of
	  parameters. This matrix includes all thinned samples.
	  }
     \item{Posterior2}{
          This is a matrix equal to \code{Posterior1}, except that rows
          are included only if stationarity (a lack of trend) is indicated 
          by the \code{Geweke.Diagnostic} for all parameters. If
	  stationarity did not occur, then this matrix is missing.
	  }
     \item{Rec.BurnIn.Thinned}{
          This is the recommended burn-in for the thinned samples, where
          the value indicates the first row that was stationary across
	  all parameters, and previous rows are discarded as burn-in.
	  Samples considered as burn-in are discarded because they do
	  not represent the target distribution and have not adequately
	  forgotten the initial value of the chain (or Markov chain, if
	  \code{Adaptive=0}).
	  }
     \item{Rec.BurnIn.UnThinned}{
          This is the recommended burn-in for all samples, in case
          thinning will not be necessary.
	  }
     \item{Rec.Thinning}{
	  This is the recommended value for the \code{Thinning} 
          argument according to the autocorrelation in the thinned
	  samples, and it is limited to the interval [1,1000].
	  }
     \item{Status}{
          This is the value in the \code{Status} argument.
          }
     \item{Summary1}{
          This is a matrix that summarizes the marginal posterior
          distributions of the parameters, deviance, and monitored
          variables over all samples in \code{Posterior1}. The following 
	  summary statistics are included: mean, standard deviation, MCSE 
          (Monte Carlo Standard Error), ESS is the effective sample size 
	  due to autocorrelation, and finally the 2.5\%, 50\%, and
	  97.5\% quantiles are reported. MCSE is essentially a standard
	  deviation around the marginal posterior mean that is due to
	  uncertainty associated with using MCMC. The acceptable size of 
          the MCSE depends on the acceptable uncertainty associated
	  around the marginal posterior mean. Laplace's Demon prefers
	  to continue updating until each MCSE is less than 6.7\% of each
	  marginal posterior standard deviation (see \code{Consort}). The 
	  The desired precision of ESS depends on the user's goal, and 
	  Laplace's Demon prefers to continue until each ESS is at least 
	  100, which should be enough to describe 95\% boundaries of an
	  approximately Gaussian distribution.
	  }
     \item{Summary2}{
          This matrix is identical to the matrix in \code{Summary1}, 
          except that it is calculated only on the stationary samples
	  found in \code{Posterior2}. If universal stationarity was not
	  estimated, then this matrix is set to missing values.
	  }
     \item{Thinned.Samples}{
          This is the number of thinned samples that were retained.
          }
     \item{Thinning}{
          This is the value of the \code{Thinning} argument.
          }
}
\references{
     Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
     Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
     Hall, London.
  
     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}

     Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). "DRAM:
     Efficient Adaptive MCMC". Statistical Computing, 16, p. 339-354.
     
     Haario, H., Saksman, E., and Tamminen, J. (2001). "An Adaptive
     Metropolis Algorithm". Bernoulli, 7, p. 223-242.

     Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". Journal of the
     American Statistical Association, 90(430), p. 773--795.

     Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
     Posterior Simulation with the Laplace-Metropolis
     Estimator". Journal of the American Statistical Association, 92,
     p. 648--655.
     
     Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller,
     E. (1953). "Equation of State Calculations by Fast Computing
     Machines". Journal of Chemical Physics, 21, p. 1087-1092.

     Mira, A. (2001). "On Metropolis-Hastings Algorithms with Delayed
     Rejection". Metron, Vol. LIX, n. 3-4, p. 231-241.
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
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("LP","tau")
parm.names <- rep(NA,J+1)
for (j in 1:J) {parm.names[j] <- paste("beta[",j,"]",sep="")}
parm.names[J+1] <- "log.tau"
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
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,tau), yhat=mu,
          parm=parm)
     return(Modelout)
     }

set.seed(666)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

##########################  Adaptive Metropolis  ##########################
Fit <- LaplacesDemon(Model, Data=MyData, Adaptive=500,
     Covar=NULL, DR=0, Initial.Values, Iterations=10000,
     Periodicity=10, Status=1000, Thinning=10)
Fit
print(Fit)
Consort(Fit)
caterpillar.plot(Fit, Parms=1:10)
BurnIn <- Fit$Rec.BurnIn.Thinned
plot(Fit, BurnIn, MyData, PDF=FALSE, Parms=Fit$Parameters)
Pred <- predict(Fit, Model, MyData)
summary(Pred)
plot(Pred)

#####################  Delayed Rejection Metropolis  ######################
Fit <- LaplacesDemon(Model, Data=MyData, Adaptive=0,
     Covar=NULL, DR=1, Initial.Values, Iterations=1000,
     Periodicity=0, Status=100, Thinning=1)

#################  Delayed Rejection Adaptive Metropolis  #################
Fit <- LaplacesDemon(Model, Data=MyData, Adaptive=500,
     Covar=NULL, DR=1, Initial.Values, Iterations=1000,
     Periodicity=100, Status=100, Thinning=1)

########################  Random-Walk Metropolis  #########################
Fit <- LaplacesDemon(Model, Data=MyData, Adaptive=0,
     Covar=NULL, DR=0, Initial.Values, Iterations=1000,
     Periodicity=0, Status=100, Thinning=1)

#End
}
\keyword{
Adaptive, Adaptive MCMC, AM, Bayesian Inference, Delayed Rejection,
DR, DRAM, DRM, Laplace Approximation, LaplacesDemon, Laplace's Demon,
Markov chain Monte Carlo, MCMC, Metropolis, Random Walk, Random-Walk,
RWM
}
