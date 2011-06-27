###########################################################################
# indmat                                                                  #
#                                                                         #
# The purpose of the indmat function is to create an indicator matrix     #
# from a vector. This function is similar to the class.ind function in    #
# the nnet package.                                                       #
###########################################################################

indmat <- function(x)
     {
     n <- length(x)
     x <- as.factor(x)
     X <- matrix(0, n, length(levels(x)) )
     X[(1:n) + n*(unclass(x)-1)] <- 1
     dimnames(X) <- list(names(x), levels(x))
     return(X)
     }

#End
