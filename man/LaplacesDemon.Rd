\name{LaplacesDemon}
\alias{LaplacesDemon}
\title{Laplace's Demon}
\description{
The \code{LaplacesDemon} function is the main function of Laplace's
Demon. Given data, a model specification, and initial values,
\code{LaplacesDemon} maximizes the logarithm of the unnormalized joint
posterior density with MCMC and provides samples of the marginal posterior
distributions, deviance, and other monitored variables.
}
\usage{
LaplacesDemon(Model, Data, Initial.Values, Covar=NULL, Iterations=100000,
     Status=1000, Thinning=100, Algorithm="RWM", Specs=NULL)
}
\arguments{
  \item{Model}{
    This required argument receives the model from a user-defined
    function. The user-defined function is where the model is
    specified. \code{LaplacesDemon} passes two arguments to the  model
    function, \code{parms} and \code{Data}, and receives five arguments
    from the model function: \code{LP} (the logarithm of the
    unnormalized joint posterior), \code{Dev} (the deviance),
    \code{Monitor} (the monitored variables), \code{yhat} (the variables
    for posterior predictive checks), and \code{parm}, the vector of
    parameters, which may be constrained in the model function.}
  \item{Data}{
    This required argument accepts a list of data. The list of data must
    contain \code{mon.names} which contains monitored variable names,
    and must contain \code{parm.names} which contains parameter
    names. The \code{\link{as.parm.names}} function may be helpful for
    preparing the data.}
  \item{Initial.Values}{
    This argument requires a vector of initial values equal in length to
    the number of parameters. Each initial value will be the starting
    point for an adaptive chain or a non-adaptive Markov chain of a
    parameter. If all initial values are set to zero, then Laplace's
    Demon will attempt to optimize the initial values with the
    \code{\link{LaplaceApproximation}} function. After Laplace's Demon
    finishes updating, it may be desired to continue updating from where
    it left off. To continue, this argument should receive the last
    iteration of the previous update. For example, if Fit is the output
    object, then \code{Initial.Values=as.initial.values(Fit)}. Initial
    values may be generated randomly with the \code{\link{GIV}}
    function.}
  \item{Covar}{
    This argument defaults to \code{NULL}, but may otherwise accept a
    \eqn{K \times K}{K x K} proposal covariance matrix for the first
    adaptation of the proposal covariances, where \eqn{K} is the number
    of dimensions (or parameters). The proposal covariance matrix will
    be re-estimated with each adaptation according to the entire history
    of all chains. When the model is updated for the first time,
    \code{Covar=NULL} should be used, unless there is a better estimate
    at the variance of each target distribution as well as the
    associated covariances. Once Laplace's Demon has finished updating,
    it may be desired to continue updating where it left off, in which
    case the proposal covariance matrix from the last run can be input
    into the next run. The covariance matrix may also be input from the
    \code{\link{LaplaceApproximation}} function, if used.}
  \item{Iterations}{
    This required argument accepts integers larger than 10, and
    determines the number of iterations that Laplace's Demon will update
    the parameters while searching for target distributions. The
    required amount of computer memory will increase with
    \code{Iterations}. If computer memory is exceeded, then all will be
    lost. The \code{\link{Combine}} function can be used later to
    combine multiple updates.}
  \item{Status}{
    This argument accepts integers between 1 and the number of
    iterations, and indicates how often the user would like the status
    of the number of iterations and proposal type (multivariate or
    single-component) printed to the screen. For example, if a model is
    updated for 1,000 iterations and \code{Status=200}, then a status
    message will be printed at the following iterations: 200, 400, 600,
    and 800.}
  \item{Thinning}{
    This argument accepts integers between 1 and the number of
    iterations, and indicates that every nth iteration will be retained,
    while the other iterations are discarded. If \code{Thinning=5}, then
    every 5th iteration will be retained. Thinning is performed to
    reduce autocorrelation and the number of marginal posterior
    samples.}
  \item{Algorithm}{
    This argument accepts the abbreviated name of the MCMC algorithm,
    which must appear in quotes. The argument defaults to RWM. Options
    include: AM, AMM, AMWG, DRAM, DRM, MWG, RAM, and RWM. Respectively,
    these algorithms are Adaptive Metropolis, Adaptive-Mixture
    Metropolis, Adaptive Metropolis-within-Gibbs, Delayed Rejection
    Adaptive Metropolis, Delayed Rejection Metropolis,
    Metropolis-within-Gibbs, Robust Adaptive Metropolis, and Random-Walk
    Metropolis.}
  \item{Specs}{
    This argument defaults to NULL, and accepts a list of
    specifications for the MCMC algorithm declared in the
    \code{Algorithm} argument. The specifications associated with each
    algorithm may be seen below in the examples. The \code{Adaptive}
    argument specifies at which iteration an adaptive algorithm begins
    to adapt, except for the AMWG algorithm, which begins to adapt
    immediately. Other adaptive algorithms depend on the sample
    covariance matrix of past iterations, and should begin later,
    depending on the size and complexity of the sample covariance
    matrix. The \code{Periodicity} argument specifies how often the
    adaptive algorithm should adapt. For example, if
    \code{Periodicity=10}, then the algorithm adapts every 10th
    iteration. A higher \code{Periodicity} is associated with an
    algorithm that runs faster, because it does not have to calculate
    adaptation as often, though the algorithm adapts less often to the
    target distributions, so it is a trade-off. It is recommended to use
    the lowest value that runs fast enough to suit the user, or provide
    sufficient adaptation. The AMM algorithm also requires \code{w}. The
    RAM algorithm also requires \code{alpha.star} [the desired
    acceptance rate(s), which can be a scalar or a vector equal in
    length to the number of targets]. Dist (which can be t for
    t-distributed or N for normally-distributed), and \code{gamma}
    (which controls the decay of adaptation). The SAMWG and SMWG
    algorithms require a T x K matrix Dyn of dynamic parameters. The
    USAMWG and USMWG algorithms require \code{Begin} to indicate the
    time-period at which to begin updating (filtering and predicting).}
}
\details{
  \code{LaplacesDemon} offers numerous MCMC algorithms for numerical
  approximation in Bayesian inference. The algorithms are

  \itemize{
    \item Adaptive Metropolis (AM)
    \item Adaptive-Mixture Metropolis (AMM)
    \item Adaptive Metropolis-within-Gibbs (AMWG)
    \item Delayed Rejection Adaptive Metropolis (DRAM)
    \item Delayed Rejection Metropolis (DRM)
    \item Metropolis-within-Gibbs (MWG)
    \item Robust Adaptive Metropolis (RAM)
    \item Random-Walk Metropolis (RWM)
    \item Sequential Adaptive Metropolis-within-Gibbs (SAMWG)
    \item Sequential Metropolis-within-Gibbs (SMWG)
    \item Updating Sequential Adaptive Metropolis-within-Gibbs (USAMWG)
    \item Updating Sequential Metropolis-within-Gibbs (USMWG)
  }

  References are provided below. The AM, DRAM, DRM, and RWM algorithms
  attempt multivariate proposals, but are modified to use single-component
  proposals if multivariate estimation is problematic. The AMM and RAM
  algorithms are currently multivariate-only, while AMWG and MWG are
  strictly single-component algorithms here. Additionally, the
  \code{\link{LaplaceApproximation}} function may be used to attempt to
  optimize initial values, which the \code{LaplacesDemon} function
  attempts automatically when all initial values are set to zero and
  sample size is at least five times the number of parameters. For more
  information, see the accompanying vignettes entitled
  "Bayesian Inference" and "LaplacesDemon Tutorial".
}
\value{
  \code{LaplacesDemon} returns an object of class \code{demonoid} that
  is a list with the following components:
  \item{Acceptance.Rate}{
    This is the acceptance rate of the MCMC algorithm, indicating
    percentage of iterations in which the proposals were accepted. The
    optimal acceptance rate varies with the number of parameters,
    ranging from 0.44 for one parameter (one IID Gaussian target
    distribution) to 0.234 for an infinite number of parameters (IID
    Gaussian target distributions), and 0.234 is approached quickly as
    the number of parameters increases.}
  \item{Adaptive}{
    This reports the value of the \code{Adaptive} argument.}
  \item{Algorithm}{
    This reports the specific algorithm used.}
  \item{Call}{
    This is the matched call of \code{LaplacesDemon}.}
  \item{Covar}{
    This stores the \eqn{d \times d}{d x d} proposal covariance matrix
    of the most recent adaptation, where \eqn{d} is the dimension or
    number of parameters or initial values. If the model is updated in
    the future, then this matrix can be used to start the next update
    where the last update left off. Only the diagonal of this matrix is
    reported in the associated \code{print} function.}
  \item{CovarDHis}{
    This \eqn{n \times d}{n x d} matrix stores the diagonal of the
    proposal covariance matrix of each adaptation in each of \eqn{n}
    rows for \eqn{d} dimensions, where the dimension is the number of
    parameters or length of the initial values vector. The proposal
    covariance matrix should change less over time.}
  \item{Deviance}{
    This is a vector of the deviance of the model, with a length equal
    to the number of thinned samples that were retained. Deviance is
    useful for considering model fit, and is equal to the sum of the
    log-likelihood for all rows in the data set, which is then
    multiplied by negative two.}
  \item{DIC1}{
    This is a vector of three values: Dbar, pD, and DIC. Dbar is the
    mean deviance, pD is a measure of model complexity indicating the
    effective number of parameters, and DIC is the Deviance Information
    Criterion, which is a model fit statistic that is the sum of Dbar
    and pD. \code{DIC1} is calculated over all retained samples. Note
    that pD is calculated as \code{var(Deviance)/2} as in Gelman et
    al. (2004).}
  \item{DIC2}{
    This is identical to \code{DIC1} above, except that it is calculated
    over only the samples that were considered by the
    \code{Geweke.Diagnostic} to be stationary for all parameters. If
    stationarity (or a lack of trend) is not estimated for all
    parameters, then \code{DIC2} is set to missing values.}
  \item{DR}{
    This reports the value of the \code{DR} argument.}
  \item{Initial.Values}{
    This is the vector of \code{Initial.Values}, which may have been
    optimized with the \code{\link{LaplaceApproximation}} function.}
  \item{Iterations}{
    This reports the number of \code{Iterations} for updating.}
  \item{LML}{
    This is an approximation of the logarithm of the marginal likelihood
    of the data (see the \code{\link{LML}} function for more
    information). \code{LML} is estimated only with stationary samples,
    and only with a non-adaptive algorithm, including Delayed Rejection
    Metropolis (DRM), Metropolis-within-Gibbs (MWG), or Random-Walk
    Metropolis (RWM). \code{LML} is estimated with nonparametric
    self-normalized importance sampling (NSIS), given LL and the
    marginal posterior samples of the parameters. \code{LML} is useful
    for comparing multiple models with the \code{\link{BayesFactor}}
    function.}
  \item{Minutes}{
    This indicates the number of minutes that \code{LaplacesDemon} was
    running, and this includes the initial checks as well as time it
    took the \code{\link{LaplaceApproximation}} function, assessing
    stationarity, effective sample size (ESS), and creating summaries.}
  \item{Model}{
    This contains the model specification \code{Model}.}
  \item{Monitor}{
    This is a vector or matrix of one or more monitored variables, which
    are variables that were specified in the \code{Model} function to be
    observed as chains (or Markov chains, if \code{Adaptive=0}), but
    that were not deviance or parameters.}
  \item{Parameters}{
    This reports the number of parameters.}
  \item{Periodicity}{
    This reports the value of the \code{Periodicity} argument.}
  \item{Posterior1}{
    This is a matrix of marginal posterior distributions composed of
    thinned samples, with a number of rows equal to the number of
    thinned samples and a number of columns equal to the number of
    parameters. This matrix includes all thinned samples.}
  \item{Posterior2}{
    This is a matrix equal to \code{Posterior1}, except that rows are
    included only if stationarity (a lack of trend) is indicated by the
    \code{\link{Geweke.Diagnostic}} for all parameters. If stationarity
    did not occur, then this matrixis missing.}
  \item{Rec.BurnIn.Thinned}{
    This is the recommended burn-in for the thinned samples, where the
    value indicates the first row that was stationary across all
    parameters, and previous rows are discarded as burn-in. Samples
    considered as burn-in are discarded because they do not represent
    the target distribution and have not adequately forgotten the
    initial value of the chain (or Markov chain, if \code{Adaptive=0}).}
  \item{Rec.BurnIn.UnThinned}{
    This is the recommended burn-in for all samples, in case thinning
    will not be necessary.}
  \item{Rec.Thinning}{
    This is the recommended value for the \code{Thinning} argument
    according to the autocorrelation in the thinned samples, and it is
    limited to the interval [1,1000].}
  \item{Status}{
    This is the value in the \code{Status} argument.}
  \item{Summary1}{
    This is a matrix that summarizes the marginal posterior
    distributions of the parameters, deviance, and monitored variables
    over all samples in \code{Posterior1}. The following summary
    statistics are included: mean, standard deviation, MCSE (Monte Carlo
    Standard Error), ESS is the effective sample size due to
    autocorrelation, and finally the 2.5\%, 50\%, and 97.5\% quantiles
    are reported. MCSE is essentially a standard deviation around the
    marginal posterior mean that is due to uncertainty associated with
    using MCMC. The acceptable size of the MCSE depends on the
    acceptable uncertainty associated around the marginal posterior
    mean. Laplace's Demon prefers to continue updating until each MCSE
    is less than 6.7\% of each marginal posterior standard deviation
    (see the \code{\link{MCSE}} and \code{\link{Consort}} functions).
    The default \code{IMPS} method is used. Next, the desired precision
    of ESS depends on the user's goal, and Laplace's Demon prefers to
    continue until each ESS is at least 100, which should be enough to
    describe 95\% boundaries of an approximately Gaussian distribution
    (see the \code{\link{ESS}} for more information).}
  \item{Summary2}{
    This matrix is identical to the matrix in \code{Summary1}, except
    that it is calculated only on the stationary samples found in
    \code{Posterior2}. If universal stationarity was not estimated, then
    this matrix is set to missing values.}
  \item{Thinned.Samples}{
    This is the number of thinned samples that were retained.}
  \item{Thinning}{
    This is the value of the \code{Thinning} argument.}
}
\references{
  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.

  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}

  Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). "DRAM:
  Efficient Adaptive MCMC". Statistical Computing, 16, p. 339-354.
     
  Haario, H., Saksman, E., and Tamminen, J. (2001). "An Adaptive
  Metropolis Algorithm". Bernoulli, 7, p. 223-242.

  Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". Journal of the
  American Statistical Association, 90(430), p. 773--795.

  Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
  Posterior Simulation with the Laplace-Metropolis Estimator". Journal
  of the American Statistical Association, 92, p. 648--655.

  Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller,
  E. (1953). "Equation of State Calculations by Fast Computing
  Machines". Journal of Chemical Physics, 21, p. 1087-1092.

  Mira, A. (2001). "On Metropolis-Hastings Algorithms with Delayed
  Rejection". Metron, Vol. LIX, n. 3-4, p. 231-241.

  Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive
  MCMC". Computational Statistics and Data Analysis, 18, p. 349--367.
  
  Rosenthal, J.S. (2007). "AMCMC: An R interface for adaptive MCMC".
  Computational Statistics and Data Analysis, 51, p. 5467--5470.

  Vihola, M. (2011). "Robust Adaptive Metropolis Algorithm with Coerced
  Acceptance Rate". Statistics and Computing. Springer, Netherlands.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{as.initial.values}},
  \code{\link{as.parm.names}},
  \code{\link{BayesFactor}},
  \code{\link{Combine}},
  \code{\link{Consort}},
  \code{\link{ESS}},
  \code{\link{Geweke.Diagnostic}},
  \code{\link{GIV}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LML}},
  \code{\link{MCSE}}.
}
\examples{
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
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
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma), yhat=mu,
          parm=parm)
     return(Modelout)
     }

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

##########################  Adaptive Metropolis  ##########################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))
Fit
print(Fit)
Consort(Fit)
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
BurnIn <- Fit$Rec.BurnIn.Thinned
plot(Fit, BurnIn, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData)
summary(Pred, Discrep="Chi-Square")
plot(Pred, Style="Covariates", Data=MyData)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="Fitted")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residuals")

######################  Adaptive-Mixture Metropolis  ######################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AMM", Specs=list(Adaptive=500, Periodicity=10, w=0.05))

###################  Adaptive Metropolis-within-Gibbs  ####################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AMWG", Specs=list(Periodicity=50))

#################  Delayed Rejection Adaptive Metropolis  #################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="DRAM", Specs=list(Adaptive=500, Periodicity=10))

#####################  Delayed Rejection Metropolis  ######################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="DRM", Specs=NULL)

#######################  Metropolis-within-Gibbs  #########################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="MWG", Specs=NULL)

######################  Robust Adaptive Metropolis  #######################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
     Periodicity=10))

########################  Random-Walk Metropolis  #########################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="RWM", Specs=NULL)

##############  Sequential Adaptive Metropolis-within-Gibbs  ##############
#NOTE: The SAMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SAMWG", Specs=list(Dyn=Dyn, Periodicity=50))

##################  Sequential Metropolis-within-Gibbs  ###################
#NOTE: The SMWG algorithm is only for state-space models (SSMs)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="SMWG", Specs=list(Dyn=Dyn))

##########  Updating Sequential Adaptive Metropolis-within-Gibbs  #########
#NOTE: The USAMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USAMWG", Specs=list(Dyn=Dyn, Periodicity=50, Fit=Fit,
#     Begin=T.m))

##############  Updating Sequential Metropolis-within-Gibbs  ##############
#NOTE: The USMWG algorithm is only for state-space model updating
#Fit <- LaplacesDemon(Model, MyData, Initial.Values, 
#     Covar=NULL, Iterations=100000, Status=100, Thinning=100,
#     Algorithm="USMWG", Specs=list(Dyn=Dyn, Fit=Fit, Begin=T.m))

#End
}
\keyword{
  Adaptive, Adaptive MCMC, AM, AMM, AMWG, Bayesian Inference,
  Delayed Rejection, DRAM, DRM, Gibbs, Laplace Approximation,
  LaplacesDemon, Laplace's Demon, Markov chain Monte Carlo, MCMC,
  Metropolis, MWG, RAM, Random Walk, Random-Walk, Robust, RWM
}
