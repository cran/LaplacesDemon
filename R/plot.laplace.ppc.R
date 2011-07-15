###########################################################################
# plot.laplace.ppc                                                        #
#                                                                         #
# The purpose of this function is to plot an object of class laplace.ppc. #
###########################################################################

plot.laplace.ppc <- function(x, Rows=NULL, PDF=FALSE,  ...)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is NULL.\n")
     if(is.null(Rows)) {Rows <- 1:nrow(x$yhat)}
     if(PDF == TRUE)
          {
          pdf("Laplace.Approximation.PPC.Plots.pdf")
          par(mfrow=c(3,3))
          }
     else {par(mfrow=c(3,3), ask=TRUE)}
     ### Plot, Posterior Predictive of yhat
     for (j in 1:length(Rows))
          {
          plot(density(x$yhat[Rows[j],]),
               main=paste("Post. Pred. Plot of yhat[", Rows[j],
                    ",]", sep=""), xlab="Value",
               sub="Black=Density, Red=y")
          polygon(density(x$yhat[Rows[j],]), col="black", border="black")
          abline(v=x$y[Rows[j]], col="red")
          }
     ### Plot Deviance
     plot(density(x$deviance), main="Deviance", xlab="Value")
     polygon(density(x$deviance), col="black", border="black")
     abline(v=0, col="red")
     ### Plot Monitors
     for (j in 1:nrow(x$monitor)) {
          plot(density(x$monitor[j,]), main=rownames(x$monitor)[j],
               xlab="Value")
          polygon(density(x$monitor[j,]), col="black", border="black")
          abline(v=0, col="red")
          }
     if(PDF == TRUE) dev.off()
     }

#End
