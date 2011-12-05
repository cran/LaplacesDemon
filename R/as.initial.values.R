###########################################################################
# as.initial.values                                                       #
#                                                                         #
# The purpose of the as.initial.values function is to retrieve the last   #
# posterior samples from an object of class demonoid or class laplace to  #
# serve as initial values for future updating.                            #
###########################################################################

as.initial.values <- function(x)
     {
     if((class(x) != "demonoid") & (class(x) != "laplace"))
          stop("x must be an object of class demonoid.")
     if(class(x) == "demonoid")
          initial.values <- as.vector(x$Posterior1[x$Thinned.Samples,])
     if(class(x) == "laplace")
          initial.values <- as.vector(x$Summary1[,"Mode"])
     return(initial.values)
     }

#End
