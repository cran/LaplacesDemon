###########################################################################
# BayesianBootstrap                                                       #
#                                                                         #
# The purpose of the BayesianBootstrap is to allow the user to sample     #
# from data for future bootstrapping.                                     #
###########################################################################

BayesianBootstrap <- function(X, n=1000, Status=NULL)
     {
     ### Initial Checks
     if(missing(X)) stop("X is a required argument.")
     if(!is.matrix(X)) X <- as.matrix(X)
     if(any(!is.finite(X))) stop("Non-finite values found in X.")
     S <- round(abs(n))
     if(S < 1) S <- 1
     if(is.null(Status)) Status <- S + 1
     else {
          Status <- round(abs(Status))
          if(Status < 1 | Status > S) Status <- S + 1}
     N <- nrow(X)
     J <- ncol(X)
     X <- as.matrix(X[order(X[,1]),])
     X.B <- matrix(NA,S,J)
     ### Bootstrap Samples
     for (s in 1:S) {
          if(s %% Status == 0) cat("\nBootstrapped Samples:", s)
          u <- c(0, sort(runif(N-1)), 1)
          g <- diff(u)
          X.B[s,] <- X[sample(1:N, 1, prob=g, replace=TRUE),]
          }
     cat("\n\nThe Bayesian Bootstrap has finished.\n\n")
     return(X.B)
     }

#End
