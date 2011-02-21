\name{BayesFactor}
\alias{BayesFactor}
\title{Bayes Factor}
\description{
This function calculates the Bayes factor of two fitted objects of class
\code{laplace} that were estimated with the \code{LaplaceApproximation}
function, and indicates the strength of evidence in favor of the
hypothesis (that the first model is better than the second model).
}
\usage{
BayesFactor(x)
}
\arguments{
     \item{x}{This is a list of two fitted objects of class \code{laplace}.}
}
\details{
Introduced by Harold Jeffreys, a 'Bayes factor' is a Bayesian
alternative to frequentist hypothesis testing that is most often used
for the comparison of multiple models by hypothesis testing, usually
to determine which model better fits the data (Jeffreys, 1961). Bayes
factors are notoriously difficult to compute, and the Bayes factor is
only defined when the marginal density of \code{y} under each model is
proper. However, the Bayes factor is easy to approximate with Laplace
Approximation (see the \code{LaplaceApproximation} function.

Hypothesis testing with Bayes factors is more robust than frequentist
hypothesis testing, since the Bayesian form avoids model selection bias,
evaluates evidence in favor the null hypothesis, includes model
uncertainty, and allows non-nested models to be compared (though of
course the model must have the same dependent variable). Also,
frequentist significance tests become biased in favor of rejecting the
null hypothesis with sufficiently large sample size.

The \code{LaplaceApproximation} function returns the \code{LML}, the
approximate logarithm of the marginal likelihood of the data, in each
fitted object of class \code{laplace}. The \code{BayesFactor} function
calculates Bayes factor \code{B} as the exponentiated difference of
\code{LML} of model 1 (M[1]) and \code{LML} of model 2 (M[2]), where
model 1 is the first fitted object in list \code{x}, and model 2 is
the second fitted object. The 'Strength of Evidence' is reported
according to the interpretation of Harold Jeffreys (Jeffreys, 1961).

The Bayes factor, \code{B}, is the posterior odds in favor of the
hypothesis divided by the prior odds in favor of the hypothesis, where
the hypothesis is usually M[1] > M[2]. For example, when \code{B=2}, the
data favor M[1] over M[2] with 2:1 odds.

The \code{BayesFactor} function currently does not work with fitted
objects of class \code{demonoid} from the \code{LaplacesDemon}
function. There are a variety of methods to approximate the logarithm of
the marginal likelihood from an object fit by MCMC. This will be
included in future releases of the \code{LaplacesDemon} package.

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
     \item{B}{This is the Bayes factor.}
     \item{Hypothesis}{
          This is the hypothesis that model 1 (M[1]) is greater than
          model 2 (M[2]).
	  }
     \item{Strength.of.Evidence}{
          This is the strength of evidence in favor of the hypothesis.
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
}
\seealso{\code{\link{LaplaceApproximation}}}
\examples{
# The following example fits a model as Fit1, then adds a predictor, and
# fits another model, Fit2. The two models are compared with a Bayes
# factor.

library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
N <- NROW(demonsnacks)
J <- 2
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,7]))
X[,2] <- CenterScale(X[,2])
mon.names <- c("tau","mu[1]")
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
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(tau,mu[1]), yhat=mu,
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
mon.names <- c("tau","mu[1]")
parm.names <- rep(NA,J+1)
for (j in 1:J) {parm.names[j] <- paste("beta[",j,"]",sep="")}
parm.names[J+1] <- "log.tau"
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