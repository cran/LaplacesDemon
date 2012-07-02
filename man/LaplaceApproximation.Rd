\name{LaplaceApproximation}
\alias{LaplaceApproximation}
\title{Laplace Approximation}
\description{
  The \code{LaplaceApproximation} function deterministically maximizes
  the logarithm of the unnormalized joint posterior density with one of
  several optimization algorithms. The goal of Laplace Approximation is
  to estimate the posterior mode and variance of each parameter. This
  function is useful for optimizing initial values and estimating a
  covariance matrix to be input into the \code{\link{LaplacesDemon}}
  function, or sometimes for model estimation in its own right.
}
\usage{
LaplaceApproximation(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="LBFGS", Samples=1000, sir=TRUE,
     Stop.Tolerance=1.0E-5)
}
\arguments{
  \item{Model}{This required argument receives the model from a
    user-defined function. The user-defined function is where the model
    is specified. \code{LaplaceApproximation} passes two arguments to
    the model function, \code{parms} and \code{Data}. For more
    information, see the \code{\link{LaplacesDemon}} function and
    ``LaplacesDemon Tutorial'' vignette.}
  \item{parm}{This argument requires a vector of initial values equal in
    length to the number of parameters. \code{LaplaceApproximation} will
    attempt to optimize these initial values for the parameters, where
    the optimized values are the posterior modes, for later use with the
    \code{\link{LaplacesDemon}} function. The \code{\link{GIV}} function
    may be used to randomly generate initial values.}
  \item{Data}{This required argument accepts a list of data. The list of
    data must include \code{mon.names} which contains monitored variable
    names, and \code{parm.names} which contains parameter
    names. \code{LaplaceApproximation} must be able to determine the
    sample size of the data, and will look for a scalar sample size
    variable \code{n} or \code{N}. If not found, it will look for
    variable \code{y} or \code{Y}, and attempt to take its number of
    rows as sample size. \code{LaplaceApproximation} needs to determine
    sample size due to the asymptotic nature of this method.}
  \item{Interval}{This argument receives an interval for estimating
    approximate gradients. The logarithm of the unnormalized joint
    posterior density of the Bayesian model is evaluated at +/- 1/2 of
    this interval. \code{Interval} defaults to 1.0E-6.}
  \item{Iterations}{This argument accepts an integer that determines the
    number of iterations that \code{LaplaceApproximation} will attempt
    to maximize the logarithm of the unnormalized joint posterior
    density. \code{Iterations} defaults to 100.
    \code{LaplaceApproximation} will stop before this number of
    iterations if the tolerance is less than or equal to the
    \code{Stop.Tolerance} criterion. The required amount of computer
    memory increases with \code{Iterations}. If computer memory is
    exceeded, then all will be lost.}
  \item{Method}{This optional argument specifies the method used for
    Laplace Approximation. The default method is \code{LBFGS}. Options
    include \code{AGA} for adaptive gradient ascent and \code{Rprop} for
    resilient backpropagation. The previous default, \code{CG}
    (conjugate gradient), has been deprecated.}
  \item{Samples}{This argument indicates the number of posterior samples
    to be taken with sampling importance resampling via the
    \code{\link{SIR}} function, which occurs only when
    \code{sir=TRUE}. Note that the number of samples should increase
    with the number and intercorrelations of the parameters.}
  \item{sir}{This logical argument indicates whether or not sampling
    importance resampling is conducted via the \code{\link{SIR}}
    function to draw independent posterior samples. This argument
    defaults to \code{TRUE}. Even when \code{TRUE}, posterior samples
    are drawn only when \code{LaplaceApproximation} has
    converged. Posterior samples are required for many other functions,
    including \code{plot.laplace} and \code{predict.laplace}. The only
    time that it is advantageous for \code{sir=FALSE} is when
    \code{LaplaceApproximation} is used to help the initial values for
    \code{LaplacesDemon}, and it is unnecessary for time to be spent on
    sampling.}
  \item{Stop.Tolerance}{This argument accepts any positive number and
    defaults to 1.0E-5. At each iteration, the square root of the sum of
    the squared differences of the logarithm of the unnormalized joint
    posterior density is calculated. If this result is less than or
    equal to the value of \code{Stop.Tolerance}, then
    \code{LaplaceApproximation} has converged to the user-specified
    tolerance and will terminate at the end of the current iteration.}
}
\details{
  The Laplace Approximation or Laplace Method is a family of asymptotic
  techniques used to approximate integrals. Laplace's method seems to
  accurately approximate unimodal posterior moments and marginal
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
  Shachter, 1994).

  The speed of Laplace Approximation depends on the optimization
  algorithm selected, and typically involves many evaluations of the
  objective function per iteration (where the AM MCMC algorithm
  evaluates once per iteration), making most of the MCMC algorithms
  faster per iteration. The attractiveness of Laplace Approximation is
  that it typically improves the objective function better than MCMC
  when the parameters are in low-probability regions (in which MCMC
  algorithms may suffer unreasonably low acceptance rates) until an
  adaptive MCMC has ``learned'' how to move better. Laplace
  Approximation is also typically faster because it is seeking
  point-estimates, rather than attempting to represent the target
  distribution with enough simulation draws. Laplace Approximation
  extends MLE, but shares similar limitations, such as its asymptotic
  nature with respect to sample size. Bernardo and Smith (2000) note
  that Laplace Approximation is an attractive numerical approximation
  algorithm, and will continue to develop.

  \code{LaplaceApproximation} seeks a global maximum of the logarithm of
  the unnormalized joint posterior density. The approach differs by
  \code{Method}. The \code{\link{LaplacesDemon}} function uses the
  \code{LaplaceApproximation} algorithm to optimize initial values,
  estimate covariance, and save time for the user.

  Most optimization algorithms assume that the logarithm of the
  unnormalized joint posterior density is defined and differentiable. An
  approximate gradient is taken for each initial value as the difference
  in the logarithm of the unnormalized joint posterior density due to a
  slight increase versus decrease in the parameter.

  When \code{Method="AGA"}, at 10 evenly-space times,
  \code{LaplaceApproximation} attempts several step sizes, which are
  also called rate parameters in other literature, and selects the best
  step size from a set of 10 fixed options. Thereafter, each iteration
  in which an improvement does not occur, the step size shrinks, being
  multiplied by 0.999.

  Gradient ascent is criticized for sometimes being relatively slow when
  close to the maximum, and its asymptotic rate of convergence is
  inferior to other methods. However, compared to other popular
  optimization algorithms such as Newton-Raphson, an advantage of the
  gradient ascent is that it works in infinite dimensions, requiring
  only sufficient computer memory. Although Newton-Raphson converges in
  fewer iterations, calculating the inverse of the negative Hessian
  matrix of second-derivatives is more computationally expensive and
  subject to singularities. Therefore, gradient ascent takes longer to
  converge, but is more generalizable.

  When \code{Method="LBFGS"}, the limited-memory BFGS
  (Broyden-Fletcher-Goldfarb-Shanno) algorithm is called in
  \code{optim}, once per iteration.

  When \code{Method="Rprop"}, the approximate gradient is taken for each
  parameter in each iteration, and its sign is compared to the
  approximate gradient in the previous iteration. A weight element in a
  weight vector is associated with each approximate gradient. A weight
  element is multiplied by 1.2 when the sign does not change, or by 0.5
  if the sign changes. The weight vector is the step size, and is
  constrained to the interval [0.001, 50], and initial weights are
  0.0125. This is the resilient backpropagation algorithm, which is
  often denoted as the ``Rprop-'' algorithm of Riedmiller (1994).

  After \code{LaplaceApproximation} finishes, due either to early
  convergence or completing the number of specified iterations, it
  approximates the Hessian matrix of second derivatives, and attempts to
  calculate the covariance matrix by taking the inverse of the negative
  of this matrix. If successful, then this covariance matrix may be
  passed to \code{\link{LaplacesDemon}}, and the diagonal of this matrix
  is the variance of the parameters. If unsuccessful, then a scaled
  identity matrix is returned, and each parameter's variance will be 1.

}
\value{
  \code{LaplaceApproximation} returns an object of class \code{laplace}
  that is a list with the following components:
  \item{Call}{This is the matched call of \code{LaplaceApproximation}.}
  \item{Converged}{This is a logical indicator of whether or not
    \code{LaplaceApproximation} converged within the specified
    \code{Iterations} according to the supplied \code{Stop.Tolerance}
    criterion. Convergence does not indicate that the global maximum has
    been found, but only that the tolerance was less than or equal to
    the \code{Stop.Tolerance} criterion.}
  \item{Covar}{This covariance matrix is the negative inverse of an
    approximate Hessian matrix. If an error is encountered in attempting
    to solve the matrix inversion, then an identity matrix is
    returned. The \code{Covar} matrix may be scaled and input into the
    \code{Covar} argument of the \code{\link{LaplacesDemon}} function
    for further MCMC estimation, or the diagonal of this matrix may be
    used to represent the posterior variance of the parameters, provided
    the algorithm converged and matrix inversion was successful. To
    scale this matrix for use with Laplace's Demon, multiply it by
    \eqn{2.38^2/d}, where \eqn{d} is the number of initial values.}
  \item{Deviance}{This is a vector of the iterative history of the
    deviance in the \code{LaplaceApproximation} function, as it sought
    convergence.}
  \item{History}{This is a matrix of the iterative history of the
    parameters in the \code{LaplaceApproximation} function, as it sought
    convergence.}
  \item{Initial.Values}{This is the vector of initial values that was
    originally given to \code{LaplaceApproximation} in the \code{parm}
    argument.}
  \item{LML}{This is an approximation of the logarithm of the marginal
    likelihood of the data (see the \code{\link{LML}} function for more
    information). When the model has converged and \code{sir=TRUE}, the
    NSIS method is used. When the model has converged and
    \code{sir=FALSE}, the LME2 (one of the Laplace-Metropolis Estimator
    functions in Laplace's Demon) method is used, which is the
    logarithmic form of equation 4 in Lewis and Raftery (1997). As a
    rough estimate of Kass and Raftery (1995), the LME-based LML is
    worrisome when the sample size of the data is less than five times
    the number of parameters, and \code{LML} should be adequate in most
    problems when the sample size of the data exceeds twenty times the
    number of parameters (p. 778). The LME is inappropriate with
    hierarchical models. However \code{LML} is estimated, it is useful
    for comparing multiple models with the \code{BayesFactor} function.}
  \item{LP.Final}{This reports the final scalar value for the logarithm
    of the unnormalized joint posterior density.}
  \item{LP.Initial}{This reports the initial scalar value for the
    logarithm of the unnormalized joint posterior density.}
  \item{Minutes}{This is the number of minutes that
    \code{LaplaceApproximation} was running, and this includes the
    initial checks as well as drawing posterior samples and creating
    summaries.}
  \item{Monitor}{ When \code{sir=TRUE}, a number of independent
    posterior samples equal to \code{Samples} is taken, and the draws
    are stored here as a matrix. The rows of the matrix are the samples,
    and the columns are the monitored variables.}
  \item{Posterior}{When \code{sir=TRUE}, a number of independent
    posterior samples equal to \code{Samples} is taken, and the draws
    are stored here as a matrix. The rows of the matrix are the samples,
    and the columns are the parameters.}
  \item{Step.Size.Final}{This is the final, scalar \code{Step.Size}
    value at the end of the \code{LaplaceApproximation} algorithm.}
  \item{Step.Size.Initial}{This is the initial, scalar \code{Step.Size}.}
  \item{Summary}{This is a summary matrix. Rows are parameters. The
    following columns are included: Mode, SD (Standard Deviation), LB
    (Lower Bound), and UB (Upper Bound). The bounds constitute a 95\% probability
    interval.}
  \item{Tolerance.Final}{This is the last \code{Tolerance} of the
    \code{LaplaceApproxiation} algorithm. It is calculated as the square
    root of the sum of the squared differences between a new and current
    vector of parameters.}
  \item{Tolerance.Stop}{This is the \code{Stop.Tolerance} criterion.}
}
\references{
  Azevedo-Filho, A. and Shachter, R. (1994). "Laplace's Method
  Approximations for Probabilistic Inference in Belief Networks with
  Continuous Variables". In "Uncertainty in Artificial Intelligence",
  Mantaras, R. and Poole, D., Morgan Kauffman, San Francisco, CA,
  p. 28--36.
     
  Bernardo, J.M. and Smith, A.F.M. (2000). "Bayesian Theory". John
  Wiley \& Sons: West Sussex, England.

  Hall, B. (2012). "Laplace's Demon", STATISTICAT, LLC.
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

  Riedmiller, M. (1994). "Advanced Supervised Learning in Multi-Layer
  Perceptrons - From Backpropagation to Adaptive Learning
  Algorithms". Computer Standards and Interfaces, 16, p. 265--278.
     
  Tierney, L. and Kadane, J.B. (1986). "Accurate Approximations for
  Posterior Moments and Marginal Densities". Journal of the American
  Statistical Association, 81(393), p. 82--86.

  Tierney, L., Kass. R., and Kadane, J.B. (1989). "Fully Exponential
  Laplace Approximations to Expectations and Variances of Nonpositive
  Functions". Journal of the American Statistical Association,
  84(407), p. 710--716.
}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{LaplacesDemon}},
  \code{\link{GIV}},
  \code{\link{LML}},
  \code{\link{optim}}, and
  \code{\link{SIR}}.
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
mon.names <- c("sigma","mu[1]")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)

############################  Initial Values  #############################
Initial.Values <- c(rep(0,J), log(1))

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
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(sigma,mu[1]), yhat=mu,
          parm=parm)
     return(Modelout)
     }

Fit <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=10000, Method="AGA")
Fit
print(Fit)
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
#plot(Fit, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData)
summary(Pred, Discrep="Chi-Square")
plot(Pred, Style="Covariates", Data=MyData)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="Fitted")
plot(Pred, Style="Jarque-Bera")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residual Density")
plot(Pred, Style="Residuals")
Levene.Test(Pred)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
#Fit$Covar is scaled (2.38^2/d) and submitted to LaplacesDemon as Covar.
#Fit$Summary[,1] is submitted to LaplacesDemon as Initial.Values.
#End
}
\keyword{Adaptive}
\keyword{Bayesian Inference}
\keyword{Gradient Ascent}
\keyword{Initial Values}
\keyword{Optimization}
\keyword{Resilient Backpropagation}