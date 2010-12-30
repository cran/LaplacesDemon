\name{LaplacesDemon}
\alias{LaplacesDemon}
\title{Laplace's Demon}
\description{
The \code{LaplacesDemon} function is the main function of Laplace's
Demon.  Given data, a model specification, and initial values,
\code{LaplacesDemon} maximizes the log of the joint posterior
distribution and provides samples of the marginal posterior
distributions, deviance, and other monitored variables.
}
\usage{
LaplacesDemon( Log.Posterior=NULL, Data=NULL, Adaptive=0, Covar=NULL,
     DR=1, Initial.Values, Iterations=1000, Periodicity=1, Status=100,
     Thinning=1)
}
\arguments{
     \item{Log.Posterior}{
          This required argument specifies the model in a user-defined
          function.
          }
     \item{Data}{
          This required argument accepts a list of data.
          } 
     \item{Adaptive}{
          This argument instructs Laplace's Demon to adapt the proposal
          covariance matrix, or proposal variance, based on the history of
          each chain, beginning with a specified iteration.  If \code{Adaptive=0}
          or is greater than the number of iterations, then adaptation will 
          not occur.  For example, if \code{Adaptive=100}, then the proposal
          variances may begin to adapt at the 100th iteration.  If the
          argument is unspecified, then it defaults to \code{Adaptive=0}, and does 
          not adapt.  Note that adaptation also depends on the \code{Periodicity}
	  argument (see below).
          }
     \item{Covar}{
          This argument defaults to \code{NULL}, but may otherwise accept a K x K
          proposal covariance matrix for the first adaptation of the proposal 
          covariances, where K is the number of dimensions (or parameters).  
          The proposal covariance matrix will be re-estimated with each
          adaptation according to the entire history of all chains.  When
          the model is updated for the first time, \code{NULL} should be used,
          unless there is a better estimate at the variance of each
          target distribution as well as the associated covariances.  Once
          Laplace's Demon has finished updating, it may be desired to
          continue updating where it left off, in which case the proposal
          covariance matrix from the last run can be input into the next run.
        }
     \item{DR}{
	  This argument indicates whether or not Delayed Rejection will
          occur, and it defaults to \code{DR=1}.  If Delayed Rejection
          does occur, then a rejected proposal will be followed by
          another proposal.  Otherwise, rejections are handled in the
          traditional way for a Metropolis sampler.
	  }
     \item{Initial.Values}{
          This argument requires a vector of initial values equal in
          length to the number of parameters.  Each initial value will be
          the starting point for an adaptive chain or a non-adaptive
	  Markov chain of a parameter.  After Laplace's Demon is finished,
          it may be desired to continue updating from where it left off.
          To continue, this argument should receive the last iteration
          of the previous update.  For example, if MyModel is the output 
          object, then 
          \code{Initial.Values=MyModel$Posterior1[NROW(MyModel$Posterior1),]}.
          }
     \item{Iterations}{
          This required argument accepts integers greater than 10, and
          determines the number of iterations that Laplace's Demon will
          update the parameters while searching for target distributions.
          }
     \item{Periodicity}{
          This argument accepts positive integers and determines the
          periodicity of adaptation in terms of iterations.  For example,
          if \code{Periodicity=10}, then adaptation will occur every 10
          iterations, and begins after the number of iterations specified
          in the \code{Adaptive} argument (above).  When \code{Periodicity=1}, it will
          adapt every iteration.  Frequent adaptation is desirable, but
          requires more time due to calculations.
          }
     \item{Status}{
          This argument accepts integers between 1 and the number of 
          iterations, and indicates how often the user would like the status
          of the number of iterations and proposal type (multivariate or 
          independent) printed to the screen.  For example, if a model is 
          updated for 1,000 iterations and \code{Status=200}, then a status message
          will be printed at the following iterations: 200, 400, 600, and
	  800.
          }
     \item{Thinning}{
          This argument accepts integers between 1 and the number of
          iterations, and indicates that every nth iteration will be
          retained, while the other iterations are discarded.  If
          \code{Thinning=5}, then every 5th iteration will be retained.  Thinning
          is performed to reduce autocorrelation and the number of
          marginal posterior samples.
	  }
}
\details{
     \code{LaplacesDemon} offers four customized MCMC algorithms for 
     numerical approximation in Bayesian inference.  The four algorithms
     are 1) Adaptive Metropolis (AM), 2) Delayed Rejection Metropolis
     (DRM), 3) Delayed Rejection Adaptive Metropolis (DRAM), and 4)
     Random-Walk Metropolis (RWM).  References are provided below.
     All algorithms attempt multivariate proposals, but are modified
     to use independent proposals if multivariate estimation is
     problematic.
}
\value{
     \code{LaplacesDemon} returns an object of class \code{demonoid}
     that is a list with the following items:
     \item{Acceptance.Rate}{
          This is the acceptance rate of the MCMC algorithm, indicating
          percentage of iterations in which the proposals were
          accepted.  The optimal acceptance rate varies with the 
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
	  This stores the proposal covariance matrix of the most recent
          adaptation.  If the model is updated in the future, then this 
          matrix can be used to start the next update where the last update
	  left off.
	  }
     \item{Deviance}{
          This is a vector of the deviance of the model, with a length
          equal to the number of thinned samples that were retained.
          Deviance is useful for considering model fit, and is equal to
	  the sum of the log-likelihood for all rows in the data set,
	  which is then multiplied by negative two.
	  }
     \item{DIC1}{
	  This is a vector of three values: Dbar, pD, and DIC.  Dbar is
          the mean deviance, pD is a measure of model complexity
	  indicating the effective number of parameters, and DIC is the
	  Deviance Information Criterion, which is a model fit statistic
	  that is the sum of Dbar and pD.  \code{DIC1} is calculated
	  over all retained samples.
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
     \item{Iterations}{
          This reports the number of \code{Iterations} for updating.
          }
     \item{Log.Posterior}{
          This contains the model specification \code{Log.Posterior}.
	  }
     \item{Minutes}{
          This indicates the number of minutes that \code{LaplacesDemon}
          was running, and this includes the initial checks as well as
	  time it took to assess burn-in, effective sample size, and
	  create summaries.
	  }
     \item{Monitor}{
          This is a vector or matrix of one or more monitored variables,
          which are variables that were specified in the
	  \code{Log.Posterior} function to be observed as chains (or
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
	  parameters.  This matrix includes all thinned samples.
	  }
     \item{Posterior2}{
          This is a matrix equal to \code{Posterior1}, except that rows
          are included only if stationarity (a lack of trend) is indicated 
          by the \code{Geweke.Diagnostic} for all parameters.  If
	  stationarity did not occur, then this matrix is missing.
	  }
     \item{Rec.BurnIn.Thinned}{
          This is the recommended burn-in for the thinned samples, where
          the value indicates the number of initial rows to discard as
	  burn-in.  Samples considered as burn-in are discarded because
	  they do not represent the target distribution and have not
          adequately forgotten the initial value of the chain (or Markov
	  chain, if \code{Adaptive=0}).
	  }
     \item{Rec.BurnIn.UnThinned}{
          This is the recommended burn-in for all samples, in case
          thinning will not be necessary.
	  }
     \item{Rec.Thinning}{
	  This is the recommended value for the \code{Thinning} 
          argument according to the autocorrelation in the thinned
	  samples, and it is limited to the interval of [1,1000].
	  }
     \item{Status}{
          This is the value in the \code{Status} argument.
          }
     \item{Summary1}{
          This is a matrix that summarizes the marginal posterior
          distributions of the parameters, deviance, and monitored
          variables over all samples in \code{Posterior1}.  The following 
	  summary statistics are included: mean, standard deviation, MCSE 
          (Monte Carlo Standard Error), Eff.Size is the effective sample
	  size due to autocorrelation, and finally the 2.5\%, 50\%, and
	  97.5\% quantiles are reported.  MCSE is essentially a standard
	  deviation around the marginal posterior mean that is due to
	  uncertainty associated with using MCMC.  The acceptable size of 
          the MCSE depends on the acceptable uncertainty associated
	  around the marginal posterior mean.  As a rough guide, it is
	  usually preferable to continue updating the model until each
	  MCSE is less than 5\% of the marginal posterior standard
	  deviation.  The desired precision of Eff.Size depends on the
	  user's goal, but usually 100 is enough to describe 95\%
	  boundaries of an approximately Gaussian distribution.
	  }
     \item{Summary2}{
          This matrix is identical to the matrix in \code{Summary1}, 
          except that it is calculated only on the stationary samples
	  found in \code{Posterior2}.  If universal stationarity was not
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
     Hall, B. (2011), "Laplace's Demon", STATISTICAT, LLC.
     URL=\url{http://www.statisticat.com/laplacesdemon.html}

     Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). DRAM:
     Efficient Adaptive MCMC. Statistical Computing, 16, p. 339-354.
     
     Haario, H., Saksman, E., and Tamminen, J. (2001). An Adaptive
     Metropolis Algorithm. Bernoulli, 7, 223-242.

     Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller,
     E. (1953). Equation of State Calculations by Fast Computing
     Machines. Journal of Chemical Physics, 21, p. 1087-1092.

     Mira, A. (2001). On Metropolis-Hastings Algorithms with Delayed
     Rejection. Metron, Vol. LIX, n. 3-4, pp. 231-241.
}

\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{\code{\link{LaplacesDemon}}}
\examples{  
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

########################  Log-Posterior Function  #########################
Log.Posterior <- function(parm, Data)
     {
     ### Prior Parameters
     beta.mu <- rep(0,J)
     beta.tau <- rep(1.0E-3,J)
     tau.alpha <- 1.0E-3
     tau.beta <- 1.0E-3
     ### Parameters
     beta <- rep(0,J)
     for (j in 1:J) {beta[j] <- parm[j]}
     tau <- exp(parm[J+1])
     ### Prior Densities
     beta.prior <- rep(0,J)
     for (j in 1:J)
          {beta.prior[j] <- dnorm(beta[j], beta.mu[j],
               1/sqrt(beta.tau[j]), log=TRUE)}
     tau.prior <- dgamma(tau, tau.alpha, tau.beta, log=TRUE)
     ### Log-Posterior
     mu <- beta %*% t(X)
     LL <- sum(dnorm(y, mu, 1/sqrt(tau), log=TRUE))
     LP <- LL + sum(beta.prior) + tau.prior
     LPout <- list(LP=LP, Dev=-2*LL, Monitor=c(tau,mu[1]), yrep=mu)
     return(LPout)
     }

set.seed(666)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

##########################  Adaptive Metropolis  ##########################
MyModel <- LaplacesDemon(Log.Posterior, Data=MyData, Adaptive=500,
     Covar=NULL, DR=0, Initial.Values, Iterations=1000,
     Periodicity=10, Status=1000, Thinning=1)
MyModel
print(MyModel)
Consort(MyModel)
#MyPred <- predict(MyModel, Log.Posterior, MyData)
#MyPred$Concordance

#####################  Delayed Rejection Metropolis  ######################
MyModel <- LaplacesDemon(Log.Posterior, Data=MyData, Adaptive=0,
     Covar=NULL, DR=1, Initial.Values, Iterations=1000,
     Periodicity=10, Status=1000, Thinning=1)
MyModel
print(MyModel)
Consort(MyModel)
#MyPred <- predict(MyModel, Log.Posterior, MyData)
#MyPred$Concordance

#################  Delayed Rejection Adaptive Metropolis  #################
MyModel <- LaplacesDemon(Log.Posterior, Data=MyData, Adaptive=500,
     Covar=NULL, DR=1, Initial.Values, Iterations=1000,
     Periodicity=10, Status=1000, Thinning=1)
MyModel
print(MyModel)
Consort(MyModel)
#MyPred <- predict(MyModel, Log.Posterior, MyData)
#MyPred$Concordance

########################  Random-Walk Metropolis  #########################
MyModel <- LaplacesDemon(Log.Posterior, Data=MyData, Adaptive=0,
     Covar=NULL, DR=0, Initial.Values, Iterations=1000,
     Periodicity=10, Status=1000, Thinning=1)
MyModel
print(MyModel)
Consort(MyModel)
#MyPred <- predict(MyModel, Log.Posterior, MyData)
#MyPred$Concordance

#################  Plot Marginal Posterior Distributions  #################
BurnIn <- MyModel$Rec.BurnIn.Thinned
plot(MyModel, BurnIn, MyData, PDF=FALSE, Parms=MyModel$Parameters)

#End
}
\keyword{
Adaptive, Adaptive MCMC, AM, Bayesian Inference, Delayed Rejection,
DR, DRAM, DRM, LaplacesDemon, Laplace's Demon, Markov chain Monte Carlo,
MCMC, Metropolis, Random Walk, Random-Walk, RWM
}
