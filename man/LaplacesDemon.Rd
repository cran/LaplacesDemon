\name{LaplacesDemon}
\alias{LaplacesDemon}
\alias{LaplacesDemon.hpc}
\title{Laplace's Demon}
\description{
  The \code{LaplacesDemon} function is the main function of Laplace's
  Demon. Given data, a model specification, and initial values,
  \code{LaplacesDemon} maximizes the logarithm of the unnormalized joint
  posterior density with MCMC and provides samples of the marginal
  posterior distributions, deviance, and other monitored variables.

  The \code{LaplacesDemon.hpc} function extends Laplace's Demon to
  parallel chains for multicore or cluster high performance computing.
}
\usage{
LaplacesDemon(Model, Data, Initial.Values, Covar=NULL, Iterations=100000,
     Status=1000, Thinning=100, Algorithm="RWM", Specs=NULL, ...)
LaplacesDemon.hpc(Model, Data, Initial.Values, Covar=NULL,
     Iterations=100000, Status=1000, Thinning=100, Algorithm="RWM",
     Specs=NULL, Chains=2, CPUs=2, Packages=NULL, Dyn.libs=NULL)
}
\arguments{
  \item{Model}{This required argument receives the model from a
    user-defined function that must be named Model. The user-defined
    function is where the model is specified. \code{LaplacesDemon}
    passes two arguments to the model function, \code{parms} and
    \code{Data}, and receives five arguments from the model function:
    \code{LP} (the logarithm of the unnormalized joint posterior),
    \code{Dev} (the deviance), \code{Monitor} (the monitored variables),
    \code{yhat} (the variables for posterior predictive checks), and
    \code{parm}, the vector of parameters, which may be constrained in
    the model function. More information on the Model specification
    function may be found in the "LaplacesDemon Tutorial" vignette, and
    the \code{\link{is.model}} function. Many examples of model
    specification functions may be found in the "Examples" vignette.}
  \item{Data}{This required argument accepts a list of data. The list of
    data must contain \code{mon.names} which contains monitored variable
    names, and must contain \code{parm.names} which contains parameter
    names. The \code{\link{as.parm.names}} function may be helpful for
    preparing the data, and the \code{\link{is.data}} function may be
    helpful for checking data.}
  \item{Initial.Values}{For \code{LaplacesDemon}, this argument requires
    a vector of initial values equal in length to the number of
    parameters. For \code{LaplacesDemon.hpc}, this argument also accepts
    a vector, in which case the same initial values will be applied to
    all parallel chains, or the argument accepts a matrix in which each
    row is a parallel chain and the number of columns is equal in length
    to the number of parameters. When a matrix is supplied for
    \code{LaplacesDemon.hpc}, each parallel chain begins with its own
    initial values that are preferably dispersed. For both
    \code{LaplacesDemon} and \code{LaplacesDemon.hpc}, each initial
    value will be the starting point for an adaptive chain or a
    non-adaptive Markov chain of a parameter. If all initial values are
    set to zero, then Laplace's Demon will attempt to optimize the
    initial values with the \code{\link{LaplaceApproximation}}
    function. After Laplace's Demon finishes updating, it may be desired
    to continue updating from where it left off. To continue, this
    argument should receive the last iteration of the previous
    update. For example, if the output object is called Fit, then
    \code{Initial.Values=as.initial.values(Fit)}. Initial values may be
    generated randomly with the \code{\link{GIV}} function.}
  \item{Covar}{This argument defaults to \code{NULL}, but may otherwise
    accept a \eqn{K \times K}{K x K} proposal covariance matrix for the
    first adaptation of the proposal covariances, where \eqn{K} is the
    number of dimensions (or parameters). The proposal covariance matrix
    will be re-estimated with each adaptation according to the entire
    history of all chains. When the model is updated for the first time,
    \code{Covar=NULL} should be used, unless there is a better estimate
    at the variance of each target distribution as well as the
    associated covariances. Once Laplace's Demon has finished updating,
    it may be desired to continue updating where it left off, in which
    case the proposal covariance matrix from the last run can be input
    into the next run. The covariance matrix may also be input from the
    \code{\link{LaplaceApproximation}} function, if used.}
  \item{Iterations}{This required argument accepts integers larger than
    10, and determines the number of iterations that Laplace's Demon
    will update the parameters while searching for target
    distributions. The required amount of computer memory will increase
    with \code{Iterations}. If computer memory is exceeded, then all
    will be lost. The \code{\link{Combine}} function can be used later
    to combine multiple updates.}
  \item{Status}{
    This argument accepts integers between 1 and the number of
    iterations, and indicates how often the user would like the status
    of the number of iterations and proposal type (for example,
    multivariate or componentwise, or mixture, or subset) printed to
    the screen. For example, if a model is updated for 1,000 iterations
    and \code{Status=200}, then a status message will be printed at the
    following iterations: 200, 400, 600, and 800 in
    \code{LaplacesDemon}. The \code{LaplacesDemon.hpc} function
    does not print the status during parallel processing.}
  \item{Thinning}{This argument accepts integers between 1 and the
    number of iterations, and indicates that every nth iteration will be
    retained, while the other iterations are discarded. If
    \code{Thinning=5}, then every 5th iteration will be
    retained. Thinning is performed to reduce autocorrelation and the
    number of marginal posterior samples.}
  \item{Algorithm}{This argument accepts the abbreviated name of the
    MCMC algorithm, which must appear in quotes. A list of MCMC
    algorithms appears below in the Details section, and the
    abbreviated name is in parenthesis.}
  \item{Specs}{This argument defaults to \code{NULL}, and accepts a list
    of specifications for the MCMC algorithm declared in the
    \code{Algorithm} argument. The specifications associated with each
    algorithm may be seen below in the examples, must appear in the
    order shown, and are described in the details section below.}
  \item{Chains}{This argument is required only for
    \code{LaplacesDemon.hpc}, and indicates the number of parallel
    chains.}
  \item{CPUs}{This argument is required only for
    \code{LaplacesDemon.hpc}, and indicates the number of central
    processing units (CPUs) of the computer or cluster. For example,
    when a user has a quad-core computer, \code{CPUs=4}.}
  \item{Packages}{This optional argument is for use only with
    \code{LaplacesDemon.hpc}, and defaults to \code{NULL}. This argument
    accepts a vector of package names to load into each parallel chain.
    If the \code{Model} specification depends on any packages, then
    these package names need to be in this vector.}
  \item{Dyn.libs}{This optional argument is for use only with
    \code{LaplacesDemon.hpc}, and defaults to \code{NULL}. This argument
    accepts a vector of the names of dynamic link libraries (shared
    objects) to load into each parallel chain. The libraries must be
    located in the working directory.}
  \item{...}{Additional arguments are unused.}
}
\details{
  \code{LaplacesDemon} offers numerous MCMC algorithms for numerical
  approximation in Bayesian inference. The algorithms are

  \itemize{
    \item Adaptive Hamiltonian Monte Carlo (AHMC)
    \item Adaptive Metropolis (AM)
    \item Adaptive Metropolis-within-Gibbs (AMWG)
    \item Adaptive-Mixture Metropolis (AMM)
    \item Componentwise Hit-And-Run Metropolis (CHARM)
    \item Delayed Rejection Adaptive Metropolis (DRAM)
    \item Delayed Rejection Metropolis (DRM)
    \item Hamiltonian Monte Carlo (HMC)
    \item Hamiltonian Monte Carlo with Dual-Averaging (HMCDA)
    \item Hit-And-Run Metropolis (HARM)
    \item Independence Metropolis (IM)
    \item Metropolis-within-Gibbs (MWG)
    \item No-U-Turn Sampler (NUTS)
    \item Random-Walk Metropolis (RWM)
    \item Reversible-Jump (RJ)
    \item Robust Adaptive Metropolis (RAM)
    \item Sequential Adaptive Metropolis-within-Gibbs (SAMWG)
    \item Sequential Metropolis-within-Gibbs (SMWG)
    \item Tempered Hamiltonian Monte Carlo (THMC)
    \item t-walk (twalk)
    \item Updating Sequential Adaptive Metropolis-within-Gibbs (USAMWG)
    \item Updating Sequential Metropolis-within-Gibbs (USMWG)
  }

  It is a goal for the documentation in the \pkg{LaplacesDemon} to be
  extensive. However, details of MCMC algorithms are best explored
  online at \url{www.bayesian-inference.com/mcmc.shtml}, as well
  as in the "LaplacesDemon Tutorial" vignette, and the "Bayesian
  Inference" vignette. Algorithm specifications (\code{Specs}) are
  listed below:

  \itemize{
    \item \code{A} is the number of initial, adaptive iterations to be
    discarded as burn-in, and is used in HMCDA and NUTS.
    \item \code{Adaptive} is the iteration in which adaptation begins,
    and is used in AM, AMM, and DRAM. These algorithms adapt according
    to an observed covariance matrix, and should sample before beginning
    to adapt.
    \item \code{alpha.star} is the desired acceptance rate(s) in RAM, and
    may be a scalar or a vector equal in length to the number of
    targets. The recommended value is \code{alpha.star=0.234}.
    \item \code{at} affects the traverse move in twalk. \code{at=6} is
    recommended. It helps when some parameters are highly correlated,
    and the correlation structure may change through the
    state-space. The traverse move is associated with an acceptance rate
    that decreases as the number of parameters increases, and is the
    reason that \code{n1} is used to select a subset of parameters each
    iteration. If adjusted, it is recommended to stay in the interval
    [2,10].
    \item \code{aw} affects the walk move in twalk, and \code{aw=1.5} is
    recommended. If adjusted, it is recommended to stay in the
    interval [0.3,2].
    \item \code{bin.n} is the scalar size parameter for a binomial prior
    distribution of model size for the RJ algorithm.
    \item \code{bin.p} is the scalar probability parameter for a
    binomial prior distribution of model size for the RJ algorithm.
    \item \code{Begin} indicates the time-period in which to begin
    updating (filtering or predicting) in the USAMWG and USMWG
    algorithms.
    \item \code{delta} is the target acceptance rate in HMCDA and
    NUTS. The recommended value is 0.65 in HMCDA and 0.6 in NUTS.
    \item \code{Dist} is the proposal distribution in RAM, and may
    either be \code{Dist="t"} for t-distributed or \code{Dist="N"} for
    normally-distributed.
    \item \code{Dyn} is a \eqn{T \times K}{T x K} matrix of dynamic
    parameters, where \eqn{T} is the number of time-periods and \eqn{K}
    is the number of dynamic parameters. \code{Dyn} is used by SAMWG,
    SMWG, USAMWG, and USMWG. Non-dynamic parameters are updated first in
    each sampler iteration, then dynamic parameters are updated in a
    random order in each time-period, and sequentially by time-period.
    \item \code{epsilon} is the step-size in AHMC, HMC, HMCDA, NUTS, and
    THMC. It is a vector equal in length to the number of parameters in
    AHMC, HMC, and THMC, and is a scalar in HMCDA and NUTS. When
    \code{epsilon=NULL} in HMCDA or NUTS (only), a reasonable initial
    value is found.
    \item \code{Fit} is an object of class \code{demonoid} in the USAMWG
    and USMWG algorithms. Posterior samples before the time-period
    specified in the \code{Begin} argument are not updated, and are used
    instead from \code{Fit}.
    \item \code{gamma} controls the decay of adaptation in RAM. It is in
    the interval (0.5,1], and 0.66 is recommended.
    \item \code{L} is a scalar number of leapfrog steps in AHMC, HMC, and
    THMC. When \code{L=1}, the algorithm reduces to Langevin Monte Carlo
    (LMC).
    \item \code{lambda} is a scalar trajectory length in HMCDA.
    \item \code{Lmax} is a scalar maximum for \code{L} (see above) in
    HMCDA.
    \item \code{mu} is a vector that is equal in length to the initial
    values. This vector will be used as the mean of the proposal
    distribution, and is usually the posterior mode of a
    previously-updated \code{\link{LaplaceApproximation}}.
    \item \code{n1} affects the size of the subset of each set of points
    to adjust, and is used in twalk. It relates to the number of
    parameters, and \code{n1=4} is recommended. If adjusted, it is
    recommended to stay in the interval [2,20].
    \item \code{parm.p} is a vector of probabilities for parameter
    selection in the RJ algorithm, and must be equal in length to
    the number of initial values.
    \item \code{Periodicity} specifies how often in iterations the
    adaptive algorithm should adapt, and is used by AHMC, AM, AMM, AMWG,
    DRAM, RAM, SAMWG, and USAMWG. If \code{Periodicity=10}, then the
    algorithm adapts every 10th iteration. A higher \code{Periodicity}
    is associated with an algorithm that runs faster, because it does
    not have to calculate adaptation as often, though the algorithm
    adapts less often to the target distributions, so it is a
    trade-off. It is recommended to use the lowest value that runs fast
    enough to suit the user, or provide sufficient adaptation.
    \item \code{selectable} is a vector of indicators of whether or not
    a parameter is selectable for variable selection in the RJ
    algorithm. Non-selectable parameters are assigned a zero, and are
    always in the model. Selectable parameters are assigned a one. This
    vector must be equal in length to the number of initial values.
    \item \code{selected} is a vector of indicators of whether or not
    each parameter is selected when the RJ algorithm begins, and
    must be equal in length to the number of initial values.
    \item \code{SIV} stands for secondary initial values and is used by
    twalk. \code{SIV} must be the same length as \code{Initial.Values},
    and each element of these two vectors must be unique from each
    other, both before and after being passed to the \code{Model}
    function. \code{SIV} defaults to \code{NULL}, in which case values
    are generated with \code{\link{GIV}}.
    \item \code{Temperature} is used in the THMC algorithm to heat up
    the momentum in the first half of the leapfrog steps, and then cool
    down the momentum in the last half. \code{Temperature} must be
    positive. When greater than 1, THMC should explore more diffuse
    distributions, and may be helpful with multimodal distributions.
    \item \code{w} is a mixture weight for the AMM algorithm. It is in
    the interval (0,1], and is recommended to use \code{w=0.05}, as per
    Roberts and Rosenthal (2009). The two mixture components are
    adaptive multivariate and static/symmetric univariate proposals. The
    mixture is determined at each iteration with mixture weight
    \code{w}. A higher value of \code{w} is associated with more
    static/symmetric univariate proposals, and a lower \code{w} is
    associated with more adaptive multivariate proposals. AMM will be
    unable to include the multivariate mixture component until it has
    accumulated some history, and models with more parameters will take
    longer to be able to use adaptive multivariate proposals.
  }
}
\value{
  \code{LaplacesDemon} returns an object of class \code{demonoid}, and
  \code{LaplacesDemon.hpc} returns an object that is a list of objects of
  class \code{demonoid.hpc}, where the number of components in the list
  is the number of parallel chains. Each object of class \code{demonoid} 
  is a list with the following components:
  \item{Acceptance.Rate}{This is the acceptance rate of the MCMC
    algorithm, indicating the percentage of iterations in which the
    proposals were accepted. The optimal acceptance rate varies with the
    number of parameters, and by algorithm. Algorithms with
    componentwise Gaussian proposals have an optimal acceptance rate of
    0.44, regardless of the number of parameters. Algorithms that update
    with multivariate Gaussian proposals tend to have an optimal
    acceptance rate that ranges from 0.44 for one parameter (one IID
    Gaussian target distribution) to 0.234 for an infinite number of
    parameters (IID Gaussian target distributions), and 0.234 is
    approached quickly as the number of parameters increases. The AHMC,
    HMC, and THMC algorithms have an optimal acceptance rate of 0.67,
    except when the algorithm specification \code{L=1}, where the
    optimal acceptance rate is 0.574.}
  \item{Adaptive}{This reports the value of the \code{Adaptive}
    argument.}
  \item{Algorithm}{This reports the specific algorithm used.}
  \item{Call}{This is the matched call of \code{LaplacesDemon}.}
  \item{Covar}{This stores the \eqn{K \times K}{K x K} proposal
    covariance matrix of the most recent adaptation, where \eqn{K} is
    the dimension or number of parameters or initial values. If the
    model is updated in the future, then this matrix can be used to
    start the next update where the last update left off. Only the
    diagonal of this matrix is reported in the associated \code{print}
    function.}
  \item{CovarDHis}{This \eqn{N \times K}{N x K} matrix stores the
    diagonal of the proposal covariance matrix of each adaptation in
    each of \eqn{N} rows for \eqn{K} dimensions, where the dimension is
    the number of parameters or length of the initial values vector. The
    proposal covariance matrix should change less over time. An
    exception is that the AHMC algorithm stores an algorithm
    specification here, which is not the diagonal of the proposal
    covariance matrix.}
  \item{Deviance}{This is a vector of the deviance of the model, with a
    length equal to the number of thinned samples that were retained.
    Deviance is useful for considering model fit, and is equal to the
    sum of the log-likelihood for all rows in the data set, which is
    then multiplied by negative two.}
  \item{DIC1}{This is a vector of three values: Dbar, pD, and DIC. Dbar
    is the mean deviance, pD is a measure of model complexity indicating
    the effective number of parameters, and DIC is the Deviance
    Information Criterion, which is a model fit statistic that is the
    sum of Dbar and pD. \code{DIC1} is calculated over all retained
    samples. Note that pD is calculated as \code{var(Deviance)/2} as in
    Gelman et al. (2004).}
  \item{DIC2}{This is identical to \code{DIC1} above, except that it is
    calculated over only the samples that were considered by the
    \code{Geweke.Diagnostic} to be stationary for all parameters. If
    stationarity (or a lack of trend) is not estimated for all
    parameters, then \code{DIC2} is set to missing values.}
  \item{DR}{
    This reports the value of the \code{DR} argument.}
  \item{Initial.Values}{This is the vector of \code{Initial.Values},
    which may have been optimized with the
    \code{\link{LaplaceApproximation}} function.}
  \item{Iterations}{This reports the number of \code{Iterations} for
    updating.}
  \item{LML}{This is an approximation of the logarithm of the marginal
    likelihood of the data (see the \code{\link{LML}} function for more
    information). \code{LML} is estimated only with stationary samples,
    and only with a non-adaptive algorithm, including Componentwise
    Hit-And-Run (CHARM), Delayed Rejection Metropolis (DRM),
    Hamiltonian Monte Carlo (HMC), Hit-And-Run (HARM),
    Metropolis-within-Gibbs (MWG), No-U-Turn Sampler (NUTS),
    Random-Walk Metropolis (RWM), Tempered Hamiltonian Monte Carlo
    (THMC) or t-walk (twalk). \code{LML} is estimated with nonparametric
    self-normalized importance sampling (NSIS), given LL and the
    marginal posterior samples of the parameters. \code{LML} is useful
    for comparing multiple models with the \code{\link{BayesFactor}}
    function.}
  \item{Minutes}{This indicates the number of minutes that
    \code{LaplacesDemon} was running, and this includes the initial
    checks as well as time it took the
    \code{\link{LaplaceApproximation}} function, assessing stationarity,
    effective sample size (ESS), and creating summaries.}
  \item{Model}{This contains the model specification \code{Model}.}
  \item{Monitor}{This is a vector or matrix of one or more monitored
    variables, which are variables that were specified in the
    \code{Model} function to be observed as chains (or Markov chains,
    if \code{Adaptive=0}), but that were not deviance or parameters.}
  \item{Parameters}{This reports the number of parameters.}
  \item{Periodicity}{This reports the value of the \code{Periodicity}
    argument.}
  \item{Posterior1}{This is a matrix of marginal posterior distributions
    composed of thinned samples, with a number of rows equal to the
    number of thinned samples and a number of columns equal to the
    number of parameters. This matrix includes all thinned samples.}
  \item{Posterior2}{This is a matrix equal to \code{Posterior1}, except
    that rows are included only if stationarity (a lack of trend) is
    indicated by the \code{\link{Geweke.Diagnostic}} for all
    parameters. If stationarity did not occur, then this matrix is
    missing.}
  \item{Rec.BurnIn.Thinned}{This is the recommended burn-in for the
    thinned samples, where the value indicates the first row that was
    stationary across all parameters, and previous rows are discarded
    as burn-in. Samples considered as burn-in are discarded because they
    do not represent the target distribution and have not adequately
    forgotten the initial value of the chain (or Markov chain, if
    \code{Adaptive=0}).}
  \item{Rec.BurnIn.UnThinned}{This is the recommended burn-in for all
    samples, in case thinning will not be necessary.}
  \item{Rec.Thinning}{This is the recommended value for the
    \code{Thinning} argument according to the autocorrelation in the
    thinned samples, and it is limited to the interval [1,1000].}
  \item{Status}{This is the value in the \code{Status} argument.}
  \item{Summary1}{This is a matrix that summarizes the marginal
    posterior distributions of the parameters, deviance, and monitored
    variables over all samples in \code{Posterior1}. The following
    summary statistics are included: mean, standard deviation, MCSE
    (Monte Carlo Standard Error), ESS is the effective sample size due
    to autocorrelation, and finally the 2.5\%, 50\%, and 97.5\%
    quantiles are reported. MCSE is essentially a standard deviation
    around the marginal posterior mean that is due to uncertainty
    associated with using MCMC. The acceptable size of the MCSE
    depends on the acceptable uncertainty associated around the
    marginal posterior mean. Laplace's Demon prefers to continue
    updating until each MCSE is less than 6.27\% of each marginal
    posterior standard deviation (see the \code{\link{MCSE}} and
    \code{\link{Consort}} functions). The default \code{IMPS} method
    is used. Next, the desired precision of ESS depends on the user's
    goal, and Laplace's Demon prefers to continue until each ESS is at
    least 100, which should be enough to describe 95\% boundaries of an
    approximately Gaussian distribution (see the \code{\link{ESS}} for
    more information).}
  \item{Summary2}{This matrix is identical to the matrix in
    \code{Summary1}, except that it is calculated only on the
    stationary samples found in \code{Posterior2}. If universal
    stationarity was not estimated for the parameters, then this matrix
    is set to missing values.}
  \item{Thinned.Samples}{This is the number of thinned samples that
    were retained.}
  \item{Thinning}{This is the value of the \code{Thinning} argument.}
}
\references{
  Christen, J.A. and Fox, C. (2010). "A General Purpose Sampling
  Algorithm for Continuous Distributions (the t-walk)". Bayesian
  Analysis, 5(2), p. 263--282.

  Duane, S., Kennedy, A.D., Pendleton, B.J., and Roweth, D. (1987).
  "Hybrid Monte Carlo". Physics Letters, B, 195, p. 216--222.
   
  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.

  Green, P.J. (1995). "Reversible Jump Markov Chain Monte Carlo
  Computation and Bayesian Model Determination". Biometrika, 82,
  p. 711--732.

  Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). "DRAM:
  Efficient Adaptive MCMC". Statistical Computing, 16, p. 339-354.
     
  Haario, H., Saksman, E., and Tamminen, J. (2001). "An Adaptive
  Metropolis Algorithm". Bernoulli, 7, p. 223-242.

  Hoffman, M.D. and Gelman. A. (2012). "The No-U-Turn Sampler: Adaptively
  Setting Path Lengths in Hamiltonian Monte Carlo". Journal of Machine
  Learning Research, p. 1--30.

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

  Smith, R.L. (1984). "Efficient Monte Carlo Procedures for Generating
  Points Uniformly Distributed Over Bounded Region". Operations
  Research, 32, p. 1296--1308.

  Vihola, M. (2011). "Robust Adaptive Metropolis Algorithm with Coerced
  Acceptance Rate". Statistics and Computing. Springer, Netherlands.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{as.initial.values}},
  \code{\link{as.parm.names}},
  \code{\link{BayesFactor}},
  \code{\link{Combine}},
  \code{\link{Consort}},
  \code{\link{ESS}},
  \code{\link{Geweke.Diagnostic}},
  \code{\link{GIV}},
  \code{\link{is.data}},
  \code{\link{is.model}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon.RAM}},
  \code{\link{LML}}, and
  \code{\link{MCSE}}.
}
\examples{
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(7,8,10)]))
J <- ncol(X)
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
PGF <- function(Data) return(c(rnormv(Data$J,0,10),
     log(rhalfcauchy(1,25))))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, y=y)

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
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

########################  Hit-And-Run Metropolis  #########################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="HARM", Specs=NULL)
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
plot(Pred, Style="ECDF")
plot(Pred, Style="Fitted")
plot(Pred, Style="Jarque-Bera")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residual Density")
plot(Pred, Style="Residuals")
Levene.Test(Pred)
Importance(Fit, Model, MyData, Discrep="Chi-Square")

##################  Adaptive Hamiltonian Monte Carlo  #####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AHMC", Specs=list(epsilon=rep(0.02, length(Initial.Values)),
#     L=2, Periodicity=10))

##########################  Adaptive Metropolis  ##########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

###################  Adaptive Metropolis-within-Gibbs  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMWG", Specs=list(Periodicity=50))

######################  Adaptive-Mixture Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="AMM", Specs=list(Adaptive=500, Periodicity=10, w=0.05))

#################  Componentwise Hit-And-Run Metropolis  ##################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="CHARM", Specs=NULL)

#################  Delayed Rejection Adaptive Metropolis  #################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRAM", Specs=list(Adaptive=500, Periodicity=10))

#####################  Delayed Rejection Metropolis  ######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="DRM", Specs=NULL)

#######################  Hamiltonian Monte Carlo  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HMC", Specs=list(epsilon=rep(0.02, length(Initial.Values)),
#     L=2))

#############  Hamiltonian Monte Carlo with Dual-Averaging  ###############
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="HMCDA", Specs=list(A=500, delta=0.65, epsilon=NULL,
#     Lmax=1000, lambda=0.1))

########################  Independence Metropolis  ########################
### Note: the mu and Covar arguments are populated from a previous Laplace
### Approximation.
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=Fit$Covar, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="IM", Specs=list(mu=Fit$Summary1[,1]))

#######################  Metropolis-within-Gibbs  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="MWG", Specs=NULL)

##########################  No-U-Turn Sampler  ############################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=100, Status=10, Thinning=1,
#     Algorithm="NUTS", Specs=list(A=50, delta=0.6, epsilon=NULL))

######################  Robust Adaptive Metropolis  #######################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#     Periodicity=10))

###########################  Reversible-Jump  #############################
#bin.n <- J-1
#bin.p <- 0.2
#parm.p <- c(1, rep(1/(J-1),(J-1)), 1)
#selectable <- c(0, rep(1,J-1), 0)
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RJ", Specs=list(bin.n=bin.n, bin.p=bin.p,
#          parm.p=parm.p, selectable=selectable,
#          selected=c(0,rep(1,J-1),0)))

########################  Random-Walk Metropolis  #########################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="RWM", Specs=NULL)

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

###################  Tempered Hamiltonian Monte Carlo  ####################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="THMC", Specs=list(epsilon=rep(0.05,length(Initial.Values)),
#     L=2, Temperature=2))

###############################  t-walk  #################################
#Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#     Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5))

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
\keyword{Adaptive}
\keyword{Adaptive Hamiltonian Monte Carlo}
\keyword{Adaptive Metropolis}
\keyword{Adaptive-Mixture Metropolis}
\keyword{Adaptive Metropolis-within-Gibbs}
\keyword{Bayesian Inference}
\keyword{Delayed Rejection Adaptive Metropolis}
\keyword{Delayed Rejection Metropolis}
\keyword{Hamiltonian Monte Carlo}
\keyword{Hit-And-Run}
\keyword{Independence Metropolis}
\keyword{MCMC}
\keyword{Metropolis-within-Gibbs}
\keyword{Multiple Chains}
\keyword{No-U-Turn Sampler}
\keyword{Optimization}
\keyword{Parallel Chains}
\keyword{Reversible-Jump}
\keyword{Robust Adaptive Metropolis}
\keyword{Random-Walk Metropolis}
\keyword{Sequential Adaptive Metropolis-within-Gibbs}
\keyword{Sequential Metropolis-within-Gibbs}
\keyword{Tempered Hamiltonian Monte Carlo}
\keyword{t-walk}
\keyword{Updating Sequential Adaptive Metropolis-within-Gibbs}
\keyword{Updating Sequential Metropolis-within-Gibbs}