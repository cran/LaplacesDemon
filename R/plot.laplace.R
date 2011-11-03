###########################################################################
# plot.laplace                                                            #
#                                                                         #
# The purpose of the plot.laplace function is to plot an object of class  #
# laplace.                                                                #
###########################################################################

plot.laplace <- function(x, Data=NULL, PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(is.null(Data)) stop("The Data argument is NULL.")
     if(is.na(x$History)) stop("There is no history to plot.")
     ### Selecting Parms
     if(is.null(Parms)) History <- x$History
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], colnames(x$History))) == 0)
               stop("Parameter in Parms does not exist.")
          keepcols <- grep(Parms[1], colnames(x$History))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i],
                         colnames(x$History))) == 0)
                         stop("Parameter in Parms does not exist.")
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
