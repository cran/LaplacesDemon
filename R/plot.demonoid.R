###########################################################################
# plot.demonoid                                                           #
#                                                                         #
# The purpose of this function is to plot an object of class demonoid.    #
###########################################################################

plot.demonoid <- function(x, BurnIn=1, Data=NULL, PDF=FALSE,
     Parms=x$Parameters, ...)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is NULL.\n")
     if(is.null(Data)) stop("The Data argument is NULL.\n")
     if(BurnIn >= NROW(x$Posterior1)) BurnIn <- 1
     if(Parms > x$Parameters) Parms <- x$Parameters
     if(PDF == TRUE)
          {
          pdf("LaplacesDemon.Plots.pdf")
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
               xlab="Value", main=Data$parm.names[j])
          polygon(density(x$Posterior1[BurnIn:x$Thinned.Samples,j]),
               col="black", border="black")
          abline(v=0, col="red", lty=2)
          z <- acf(x$Posterior1[BurnIn:x$Thinned.Samples,j], plot=FALSE)
          se <- 1/sqrt(length(x$Posterior1[BurnIn:x$Thinned.Samples,j]))
          plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
               main=Data$parm.names[j], xlab="Lag", ylab="Correlation")
          abline(h=(2*se), col="red", lty=2)
          abline(h=(-2*se), col="red", lty=2)
          }
     ### Plot Deviance
     plot(BurnIn:length(x$Deviance),
          x$Deviance[BurnIn:length(x$Deviance)],
          type="l", xlab="Iterations", ylab="Value", main="Deviance")
     panel.smooth(BurnIn:length(x$Deviance),
          x$Deviance[BurnIn:length(x$Deviance)], pch="")
     plot(density(x$Deviance[BurnIn:length(x$Deviance)]),
          xlab="Value", main="Deviance")
     polygon(density(x$Deviance[BurnIn:length(x$Deviance)]), col="black",
               border="black")
     abline(v=0, col="red", lty=2)
     z <- acf(x$Deviance[BurnIn:length(x$Deviance)], plot=FALSE)
     se <- 1/sqrt(length(x$Deviance[BurnIn:length(x$Deviance)]))
     plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
          main="Deviance", xlab="Lag", ylab="Correlation")
     abline(h=(2*se), col="red", lty=2)
     abline(h=(-2*se), col="red", lty=2)
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
               main=Data$mon.names[j])
          panel.smooth(BurnIn:nn, x$Monitor[BurnIn:nn,j], pch="")
          plot(density(x$Monitor[BurnIn:nn,j]),
               xlab="Value", main=Data$mon.names[j])
          polygon(density(x$Monitor[BurnIn:nn,j]), col="black",
               border="black")
          abline(v=0, col="red", lty=2)
          z <- acf(x$Monitor[BurnIn:nn,j], plot=FALSE)
          se <- 1/sqrt(length(x$Monitor[BurnIn:nn,j]))
          plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
               main=Data$mon.names[j], xlab="Lag", ylab="Correlation")
          abline(h=(2*se), col="red", lty=2)
          abline(h=(-2*se), col="red", lty=2)
          }
     if(PDF == TRUE) dev.off()
     }

#End
