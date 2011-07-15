###########################################################################
# plot.demonoid.ppc                                                       #
#                                                                         #
# The purpose of this function is to plot an object of class demonoid.ppc.#
###########################################################################

plot.demonoid.ppc <- function(x, Rows=NULL, PDF=FALSE,  ...)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is NULL.\n")
     if(is.null(Rows)) Rows <- 1:nrow(x$yhat)
     if(PDF == TRUE)
          {
          pdf("Laplace.Demon.PPC.Plots.pdf")
          par(mfrow=c(3,3))
          }
     else {par(mfrow=c(3,3), ask=TRUE)}
     ### Plot
     for (j in 1:length(Rows))
          {
          plot(density(x$yhat[Rows[j],]),
               main=paste("Post. Pred. Plot of yhat[", Rows[j],
                    ",]", sep=""), xlab="Value",
               sub="Black=Density, Red=y")
          polygon(density(x$yhat[Rows[j],]), col="black", border="black")
          abline(v=x$y[Rows[j]], col="red")
          }
     if(PDF == TRUE) dev.off()
     }

#End
