###########################################################################
# plot.demonoid                                                           #
#                                                                         #
# The purpose of this function is to plot an object of class demonoid.    #
###########################################################################

plot.demonoid <- function(x, BurnIn=1, Data=NULL, PDF=FALSE,
     Parms=NULL, ...)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is NULL.\n")
     if(is.null(Data)) stop("The Data argument is NULL.\n")
     if(BurnIn >= NROW(x$Posterior1)) BurnIn <- 1
     if(is.null(Parms)) Posterior <- x$Posterior1
     if(!is.null(Parms)) {
          for (i in 1:length(Parms)) {
               if(i == 1) {keepcols <- grep(Parms[i], fixed=TRUE,
                    colnames(x$Posterior1))}
               if(i > 1) {
                    newcols <- grep(Parms[i], fixed=TRUE,
                         colnames(x$Posterior1))
                    keepcols <- c(keepcols, newcols)
                    }
               }
          Posterior <- as.matrix(x$Posterior1[,keepcols])
          colnames(Posterior) <- colnames(x$Posterior1)[keepcols]
          }
     if(PDF == TRUE)
          {
          pdf("LaplacesDemon.Plots.pdf")
          par(mfrow=c(3,3))
          }
     if(PDF == FALSE) par(mfrow=c(3,3), ask=TRUE)
     ### Plot Parameters
     for (j in 1:NCOL(Posterior))
          {
          plot(BurnIn:x$Thinned.Samples,
               Posterior[BurnIn:x$Thinned.Samples,j],
               type="l", xlab="Iterations", ylab="Value",
               main=Data$parm.names[j])
          panel.smooth(BurnIn:x$Thinned.Samples,
               Posterior[BurnIn:x$Thinned.Samples,j], pch="")
          plot(density(Posterior[BurnIn:x$Thinned.Samples,j]),
               xlab="Value", main=Data$parm.names[j])
          polygon(density(Posterior[BurnIn:x$Thinned.Samples,j]),
               col="black", border="black")
          abline(v=0, col="red", lty=2)
          ### Only plot an ACF if there's > 1 unique values
          if(length(unique(Posterior[BurnIn:x$Thinned.Samples,j])) > 1) {
               z <- acf(Posterior[BurnIn:x$Thinned.Samples,j], plot=FALSE)
               se <- 1/sqrt(length(Posterior[BurnIn:x$Thinned.Samples,j]))
               plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
                    main=Data$parm.names[j], xlab="Lag", ylab="Correlation")
               abline(h=(2*se), col="red", lty=2)
               abline(h=(-2*se), col="red", lty=2)
               }
          if(length(unique(Posterior[BurnIn:x$Thinned.Samples,j])) == 1) {
               plot(0,0, main=paste(Data$parm.names[j], "is a constant."))
               }
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
          ### Only plot an ACF if there's > 1 unique values
          if(length(unique(x$Monitor[BurnIn:nn,j])) > 1) {
               z <- acf(x$Monitor[BurnIn:nn,j], plot=FALSE)
               se <- 1/sqrt(length(x$Monitor[BurnIn:nn,j]))
               plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
                    main=Data$mon.names[j], xlab="Lag", ylab="Correlation")
               abline(h=(2*se), col="red", lty=2)
               abline(h=(-2*se), col="red", lty=2)
               }
          if(length(unique(x$Monitor[BurnIn:nn,j])) == 1) {
               plot(0,0, main=paste(Data$mon.names[j], "is a constant."))
               }
          }
     ### Proposal Variance (Adaptive Algorithms only)
     if(nrow(x$CovarDHis) > 1) {
          diff <- x$CovarDHis[-1,]
          adaptchange <- matrix(NA, nrow(diff), 3)
          for (i in 2:nrow(x$CovarDHis)) {
               diff[i-1,] <- abs(x$CovarDHis[i,] - x$CovarDHis[i-1,])}
          for (i in 1:nrow(diff)) {
               adaptchange[i,1] <- quantile(diff[i,], probs=0.025)
               adaptchange[i,2] <- quantile(diff[i,], probs=0.500)
               adaptchange[i,3] <- quantile(diff[i,], probs=0.975)}
          plot(adaptchange[,2], ylim=c(min(adaptchange), max(adaptchange)),
               type="l", col="red", xlab="Adaptations",
               ylab="Absolute Difference", main="Proposal Variance",
               sub="Median=Red, 95% Bounds=Gray")
          lines(adaptchange[,1], col="gray")
          lines(adaptchange[,3], col="gray")
          }
     if(PDF == TRUE) dev.off()
     }

#End
