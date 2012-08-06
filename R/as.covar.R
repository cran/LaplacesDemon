###########################################################################
# as.covar                                                                #
#                                                                         #
# The purpose of the as.covar function is to retrieve the covariance      #
# matrix from an object of class demonoid, demonoid.hpc, laplace, or pmc, #
# or in the case of an object of class pmc with mixture components, to    #
# retrieve multiple covariance matrices.                                  #
###########################################################################

as.covar <- function(x)
     {
     if(!identical(class(x), "demonoid") &
        !identical(class(x), "demonoid.hpc") &
        !identical(class(x), "laplace") &
        !identical(class(x), "pmc"))
          stop("The class of x is unknown.")
     if(identical(class(x), "demonoid"))
          covar <- x$Covar
     else if(identical(class(x), "demonoid.hpc")) {
          Chains <- length(x)
          Deviance <- list()
          for (i in 1:Chains) {Deviance[[i]] <- x[[i]][["Deviance"]]}
          j <- which.min(sapply(Deviance, function(x)
               {min(x[length(x)])}))
          cat("\nChain",j,"has the lowest deviance.\n")
          covar <- x[[j]]$Covar}
     else if(identical(class(x), "laplace"))
          covar <- x$Covar
     else covar <- x$Covar[,,x$Iterations,]
     return(covar)
     }

#End
