\name{BayesFactor}
\alias{BayesFactor}
\title{Bayes Factor}
\description{
This function calculates Bayes factors for two or more fitted objects
of class \code{demonoid} or  \code{laplace} that were estimated
respectively with the \code{LaplacesDemon} or \code{LaplaceApproximation}
functions, and indicates the strength of evidence in favor of the
hypothesis (that each model is better than another model).
}
\usage{
BayesFactor(x)
}
\arguments{
    \item{x}{This is a list of two or more fitted objects of class
    \code{demonoid} or \code{laplace}.}
}
\details{
Introduced by Harold Jeffreys, a 'Bayes factor' is a Bayesian
alternative to frequentist hypothesis testing that is most often used
for the comparison of multiple models by hypothesis testing, usually
to determine which model better fits the data (Jeffreys, 1961). Bayes
factors are notoriously difficult to compute, and the Bayes factor is
only defined when the marginal density of \code{y} under each model is
proper. However, the Bayes factor is easy to approximate with the
Laplace-Metropolis estimator (Lewis and Raftery, 1997), and is
estimated this way in both the \code{LaplaceApproximation} and
\code{LaplacesDemon} functions.

Hypothesis testing with Bayes factors is more robust than frequentist
hypothesis testing, since the Bayesian form avoids model selection bias,
evaluates evidence in favor of the null hypothesis, includes model
uncertainty, and allows non-nested models to be compared (though of
course the model must have the same dependent variable). Also,
frequentist significance tests become biased in favor of rejecting the
null hypothesis with sufficiently large sample size.

The \code{LaplaceApproximation} and \code{LaplacesDemon} functions
each return the \code{LML}, the approximate logarithm of the marginal
likelihood of the data, in each fitted object of class \code{laplace}
and \code{demonoid}. The \code{BayesFactor} function calculates matrix
\code{B}, a matrix of Bayes factors, where each element of matrix
\code{B} is a comparison of two models. Each Bayes factor is calculated
as the exponentiated difference of \code{LML} of model 1 (M[1]) and
\code{LML} of model 2 (M[2]), and the hypothesis for each element of
matrix \code{B} is that the model associated with the row is greater
than the model associated with the column. For example, element
\code{B[3,2]} is the Bayes factor that model 3 is greater than model 2.
The 'Strength of Evidence' aids in the interpretation (Jeffreys, 1961).

Each Bayes factor, \code{B}, is the posterior odds in favor of the
hypothesis divided by the prior odds in favor of the hypothesis, where
the hypothesis is usually M[1] > M[2]. For example, when
\code{B[3,2]=2}, the data favor M[3] over M[2] with 2:1 odds.

Gelman finds Bayes factors generally to be irrelevant, because they
compute the relative probabilities of the models conditional on one
of them being true. Gelman prefers approaches that measure the
distance of the data to each of the approximate models (Gelman et al.,
2004, p. 180), such as with posterior predictive checks (see the
\code{predict.laplace} function in the context of Laplace Approximation,
or the \code{predict.demon} function in the context of MCMC). Kass et
al. (1995) asserts this can be done without assuming one model is the
true model.
}
\value{
     \code{BayesFactor} returns a list with the following components:
     \item{B}{This is a matrix of Bayes factors.}
     \item{Hypothesis}{
          This is the hypothesis, and is stated as 'row > column',
	  indicating that the model associated with the row of an
	  element in matrix \code{B} is greater than the model
	  associated with the column of that element.
	  }
     \item{Strength.of.Evidence}{
          This is the strength of evidence in favor of the hypothesis.
          }
     \item{Posterior.Probability}{
          This is a vector of the posterior probability of each model,
          given flat priors.
	  }
}
\references{
     Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
     Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
     Hall, London.
  
     Jeffreys, H. (1961). "Theory of Probability, Third Edition". Oxford
     University Press: Oxford, England.
  
     Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". Journal of the
     American Statistical Association, 90(430), p. 773--795.

     Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
     Posterior Simulation with the Laplace-Metropolis
     Estimator". Journal of the American Statistical Association, 92,
     p. 648--655.
}
\seealso{\code{\link{LaplaceApproximation}}}
\examples{
# The following example fits a model as Fit1, then adds a predictor, and
# fits another model, Fit2. The two models are compared with Bayes
# factors.

library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
N <- NROW(demonsnacks)
J <- 2
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,7]))
X[,2] <- CenterScale(X[,2])
mon.names <- c("LP","sigma")
parm.names <- parm.names(list(beta=rep(0,J), log.sigma=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)

############################  Initial Values  #############################
Initial.Values <- c(rep(0,J), log(1))

##########################  Model Specification  ##########################
Model <- function(parm, Data)
     {
     ### Parameters
     beta <- parm[1:Data$J]
     sigma <- exp(parm[Data$J+1])
     ### Log(Prior Densities)
     beta.prior <- sum(dnorm(beta, 0, sqrt(1000), log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     ### Log-Posterior
     mu <- tcrossprod(beta, Data$X)
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma), yhat=mu,
          parm=parm)
     return(Modelout)
     }

########################  Laplace Approximation  ##########################
Fit1 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=1000)
Fit1

##############################  Demon Data  ###############################
data(demonsnacks)
N <- NROW(demonsnacks)
J <- 3
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(7,8)]))
X[,2] <- CenterScale(X[,2])
X[,3] <- CenterScale(X[,3])
mon.names <- c("sigma","mu[1]")
parm.names <- parm.names(list(beta=rep(0,J), log.sigma=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)

############################  Initial Values  #############################
Initial.Values <- c(rep(0,J), log(1))

########################  Laplace Approximation  ##########################
Fit2 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=1000)
Fit2

#############################  Bayes Factor  ##############################
Model.list <- list(M1=Fit1, M2=Fit2)
BayesFactor(Model.list)
}
\keyword{Bayes factor, hypothesis testing, model selection}