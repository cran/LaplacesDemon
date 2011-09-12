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
     if(is.na(x$History)) stop("There is no history to plot.\n")
     ### Selecting Parms
     if(is.null(Parms)) History <- x$History
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          keepcols <- grep(Parms[1], colnames(x$History))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    keepcols <- c(keepcols, grep(Parms[i],
                         colnames(x$History)))}}
          History <- as.matrix(x$History[,keepcols])
          colnames(History) <- colnames(x$History)[keepcols]
          }
     if(PDF == TRUE)
          {
          pdf("LaplaceApproximation.Plots.pdf")
          par(mfrow=c(3,3))
          }
     else {par(mfrow=c(3,3), ask=TRUE)}
     ### Plot Parameters
     for (j in 1:NCOL(History))
          {
          plot(1:NROW(History), History[,j],
               type="l", xlab="Iterations", ylab="Value",
               main=colnames(History)[j])
          }
     ### Plot Deviance
     plot(1:length(x$Deviance), x$Deviance, type="l", xlab="Iterations",
          ylab="Value", main="Deviance")
     if(PDF == TRUE) dev.off()
     }

#End
