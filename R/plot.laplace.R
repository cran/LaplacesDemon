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
     if(is.null(Parms)) History <- x$History
     if(!is.null(Parms)) {
          for (i in 1:length(Parms)) {
               if(i == 1) {keepcols <- grep(Parms[i], fixed=TRUE,
                    colnames(x$History))}
               if(i > 1) {
                    newcols <- grep(Parms[i], fixed=TRUE,
                         colnames(x$History))
                    keepcols <- c(keepcols, newcols)
                    }
               }
          History <- as.matrix(x$History[,keepcols])
          colnames(History) <- colnames(x$History)[keepcols]
          }
     if(PDF == TRUE)
          {
          pdf("LaplaceApproximation.Plots.pdf")
          par(mfrow=c(3,3))
          }
     if(PDF == FALSE) par(mfrow=c(3,3), ask=TRUE)
     ### Plot Parameters
     for (j in 1:NCOL(History))
          {
          plot(1:NROW(History), History[,j],
               type="l", xlab="Iterations", ylab="Value",
               main=Data$parm.names[j])
          }
     ### Plot Deviance
     plot(1:length(x$Deviance), x$Deviance, type="l", xlab="Iterations",
          ylab="Value", main="Deviance")
     if(PDF == TRUE) dev.off()
     }

#End
