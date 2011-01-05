###########################################################################
# plot.demonoid                                                           #
#                                                                         #
# The purpose of this function is to plot an object of class demonoid.ppc.#
###########################################################################

plot.demonoid.ppc <- function(x, Rows=NULL, PDF=FALSE,  ...)
     {
     ### Initial Checks
     if(is.null(x)) {cat("ERROR: The x argument is empty.\n")}
     if(is.null(Rows)) {Rows <- 1:nrow(x$yhat)}
     if(PDF == TRUE)
          {
          pdf("Laplace.Demon.PPC.Plots.pdf")
          par(mfrow=c(3,3))
          }
     if(PDF == FALSE) par(mfrow=c(3,3), ask=TRUE)
     ### Plot
     for (j in 1:length(Rows))
          {
          plot(density(x$yhat[Rows[j],]),
               main=paste("Post. Pred. Plot of yhat[", Rows[j],
                    ",]", sep=""), sub="Black=Density, Blue=y")
          abline(v=x$y[Rows[j]], col="blue")
          }
     if(PDF == TRUE) dev.off()
     }
