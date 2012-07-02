###########################################################################
# as.initial.values                                                       #
#                                                                         #
# The purpose of the as.initial.values function is to retrieve the last   #
# posterior samples from an object of class demonoid, demonoid.hpc, or    #
# laplace to serve as initial values for future updating.                 #
###########################################################################

as.initial.values <- function(x)
     {
     if(!identical(class(x), "demonoid") &
        !identical(class(x), "demonoid.hpc") &
        !identical(class(x), "laplace"))
          stop("The class of x is unknown.")
     if(identical(class(x), "demonoid"))
          initial.values <- as.vector(x$Posterior1[x$Thinned.Samples,])
     if(identical(class(x), "demonoid.hpc")) {
          Chains <- length(x)
          Deviance <- list()
          for (i in 1:Chains) {Deviance[[i]] <- x[[i]][["Deviance"]]}
          j <- which.min(sapply(Deviance, function(x)
               {min(x[length(x)])}))
          cat("\nChain",j,"has the lowest deviance.\n")
          initial.values <- as.vector(x[[j]][["Posterior1"]][x[[j]][["Thinned.Samples"]],])}
     if(identical(class(x), "laplace"))
          initial.values <- as.vector(x$Summary1[,"Mode"])
     return(initial.values)
     }

#End
