###########################################################################
# plot.demonoid                                                           #
#                                                                         #
# The purpose of this function is to plot an object of class demonoid.    #
###########################################################################

plot.demonoid <- function(x, BurnIn=1, Data=NULL, PDF=FALSE,
     Parms=x$Parameters, ...)
     {
     ### Initial Checks
     if(is.null(Data)) {cat("ERROR: The Data argument is empty.\n")}
     if(BurnIn >= NROW(x$Posterior1)) BurnIn <- 1
     if(Parms > x$Parameters) Parms <- x$Parameters
     if(PDF == TRUE)
          {
          pdf("Laplace.Demon.Plots.pdf")
          par(mfrow=c(3,3))
          }
     if(PDF == FALSE) par(mfrow=c(3,3), ask=TRUE)
     ### Plot Parameters
     for (j in 1:Parms)
          {
          plot(BurnIn:x$Thinned.Samples,
               x$Posterior1[BurnIn:x$Thinned.Samples,j],
               type="l", xlab="Iterations", ylab="Value",
               main=Data$parm.names[j])
          panel.smooth(BurnIn:x$Thinned.Samples,
               x$Posterior1[BurnIn:x$Thinned.Samples,j], pch="")
          plot(density(x$Posterior1[BurnIn:x$Thinned.Samples,j]),
               main=Data$parm.names[j])
          acf(x$Posterior1[BurnIn:x$Thinned.Samples,j],
               main=Data$parm.names[j])
          }
     ### Plot Deviance
     plot(BurnIn:length(x$Deviance),
          x$Deviance[BurnIn:length(x$Deviance)],
          type="l", xlab="Iterations", ylab="Value", main="Deviance")
     panel.smooth(BurnIn:length(x$Deviance),
          x$Deviance[BurnIn:length(x$Deviance)], pch="")
     plot(density(x$Deviance[BurnIn:length(x$Deviance)]),
          main="Deviance")
     acf(x$Deviance[BurnIn:length(x$Deviance)],
          main="Deviance")
     ### Plot Monitored Variables
     if (is.vector(x$Monitor))
          {J <- 1; nn <- length(x$Monitor)}
     if (is.matrix(x$Monitor))
          {J <- NCOL(x$Monitor); nn <- NROW(x$Monitor)}
     for (j in 1:J)
          {
          plot(BurnIn:nn,
               x$Monitor[BurnIn:nn,j],
               type="l", xlab="Iterations", ylab="Value",
               main="Monitor")
          panel.smooth(BurnIn:nn, x$Monitor[BurnIn:nn,j], pch="")
          plot(density(x$Monitor[BurnIn:nn,j]),
               main="Monitor")
          acf(x$Monitor[BurnIn:nn,j],
               main="Monitor")
          }
     if(PDF == TRUE) dev.off()
     }

#End
