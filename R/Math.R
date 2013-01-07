###########################################################################
# Math                                                                    #
#                                                                         #
# This is a collection of functions to facilitate math.                   #
###########################################################################

logadd <- function (x, add=TRUE) 
     {
     x <- as.vector(x)
     x <- sort(x[is.finite(x)], decreasing=TRUE)
     x <- c(x[1], x[which(x != x[1])])
     if(length(x) == 1) return(x)
     n <- length(x)
     if(add == TRUE)
          z <- x[1] + log(1 + sum(exp(x[-1] - x[1])))
     else 
          z <- x[1] + sum(log(1 - exp(x[-1] - x[1])))
     return(z)
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
