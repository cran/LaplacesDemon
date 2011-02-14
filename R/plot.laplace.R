###########################################################################
# plot.laplace                                                            #
#                                                                         #
# The purpose of this function is to plot an object of class laplace.     #
###########################################################################

plot.laplace <- function(x, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is NULL.\n")
     if(is.null(Data)) stop("The Data argument is NULL.\n")
     if(is.null(Parms)) Parms <- length(x$Initial.Values)
     if(Parms > length(x$Initial.Values)) Parms <- length(x$Initial.Values)
     if(PDF == TRUE)
          {
          pdf("LaplaceApproximation.Plots.pdf")
          par(mfrow=c(3,3))
          }
     if(PDF == FALSE) par(mfrow=c(3,3), ask=TRUE)
     ### Plot Parameters
     for (j in 1:Parms)
          {
          plot(1:NROW(x$History), x$History[,j],
               type="l", xlab="Iterations", ylab="Value",
               main=Data$parm.names[j])
          }
     ### Plot Deviance
     plot(1:length(x$Deviance), x$Deviance, type="l", xlab="Iterations",
          ylab="Value", main="Deviance")
     if(PDF == TRUE) dev.off()
     }

#End
