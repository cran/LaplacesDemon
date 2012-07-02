###########################################################################
# Math                                                                    #
#                                                                         #
# This is a collection of functions to facilitate math.                   #
###########################################################################

logadd <- function(x, y, add=TRUE)
     {
     if(!is.vector(x)) x <- as.vector(x)
     if(!is.vector(y)) y <- as.vector(y)
     if(add == TRUE) out <- log(x) - log(1 + exp(log(y) + log(x)))
     else out <- log(x) + log(1 - exp(log(y) - log(x)))
     return(out)
     }
partial <- function(Model, parm, Data, Interval=1e-6)
     {
     parm.len <- length(parm)
     parm.int <- {diag(parm.len) * Interval} / 2
     parm.int[which(!is.finite(parm.int))] <- 0
     out <- rep(0, parm.len)
     for (i in 1:parm.len)
          {
          high <- Model(parm + parm.int[,i], Data)[[1]]
          low <- Model(parm - parm.int[,i], Data)[[1]]
          out[i] <- {high - low} / Interval
          if(!is.finite(high) | !is.finite(low)) out[i] <- 0
          }
     return(out)
     }

#End
