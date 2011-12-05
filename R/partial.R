###########################################################################
# partial                                                                 #
#                                                                         #
# The purpose of the partial function is to estimate a vector of first    #
# partial derivatives with finite differencing.                           #
###########################################################################

partial <- function(Model, parm, Data, Interval=1e-6)
     {
     parm.len <- length(parm)
     parm.int <- {diag(parm.len) * Interval} / 2
     parm.int <- ifelse(!is.finite(parm.int), 0, parm.int)
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
