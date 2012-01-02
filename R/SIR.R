###########################################################################
# Sampling Importance Resampling (SIR)                                    #
#                                                                         #
# The purpose of the SIR function is to perform sampling importance       #
# re-sampling, usually to draw samples from the posterior as output from  #
# LaplaceApproxmation function. This function is similar to the sir       #
# function in the LearnBayes package.                                     #
###########################################################################

SIR <- function(Model, Data, mu, Sigma, n=1000) 
     {
     if(missing(Model)) stop("The Model function is required.")
     if(missing(Data)) stop("The Data function is required.")
     if(missing(mu)) stop("The mu function is required.")
     if(!is.vector(mu)) mu <- as.vector(mu)
     if(missing(Sigma)) stop("The Sigma function is required.")
     if(!is.symmetric.matrix(Sigma)) Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma)) Sigma <- as.positive.definite(Sigma)
     if(length(mu) != nrow(Sigma)) stop("mu and Sigma are incompatible.")
     ### Sampling
     k <- length(mu)
     theta <- rmvn(n, mu, Sigma)
     theta <- ifelse(!is.finite(theta), 0, theta)
     colnames(theta) <- Data$parm.names
     ### Importance
     lf <- matrix(0, c(dim(theta)[1], 1))
     for (i in 1:n) {
          mod <- Model(theta[i,], Data)
          lf[i] <- mod[[1]]
          theta[i,] <- mod[[5]]
          }
     lp <- dmvn(theta, mu, Sigma, log=TRUE)
     md <- max(lf - lp)
     w <- exp(lf - lp - md)
     if(any(is.infinite(w))) {
          w <- ifelse(w == Inf, max(w[is.finite(w)]), w)
          w <- ifelse(w == -Inf, min(w[is.finite(w)]), w)}
     probs <- w / sum(w)
     ### Resampling
     options(warn=-1)
     test <- try(indices <- sample(1:n, size=n, replace=TRUE,
          prob=probs), silent=TRUE)
     options(warn=0)
     if(class(test) == "try-error") indices <- 1:n
     if(k > 1) theta <- theta[indices,]
     else theta <- theta[indices]
     return(theta)
     }

#End
