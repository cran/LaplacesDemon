###########################################################################
# predict.laplace                                                         #
#                                                                         #
# The purpose of the predict.laplace function is to predict y[new] or     #
# y[rep], and later provide posterior predictive checks for objects of    #
# class laplace.                                                          #
###########################################################################

predict.laplace <- function(object, Model, Data, ...)
     {
     ### Initial Checks
     if(missing(object)) stop("The object argument is required.")
     if(object$Converged == FALSE)
          stop("LaplaceApproximation did not converge.")
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data$y) & is.null(Data$Y)) stop("Data must have y or Y.")
     if(!is.null(Data$y)) y <- as.vector(Data$y)
     if(!is.null(Data$Y)) y <- as.vector(Data$Y)
     ### p(y[rep] | y), Deviance, and Monitors
     Dev <- rep(NA, nrow(object$Posterior))
     monitor <- matrix(NA, length(Data$mon.names), nrow(object$Posterior))
     lengthcomp <- as.vector(Model(object$Posterior[1,], Data)[[4]])
     if(!identical(length(lengthcomp), length(y)))
          stop("y and yhat differ in length.")
     yhat <- matrix(NA, length(y), nrow(object$Posterior))
     for (i in 1:nrow(object$Posterior)) {
          mod <- Model(object$Posterior[i,], Data)
          Dev[i] <- as.vector(mod[[2]])
          monitor[,i] <- as.vector(mod[[3]])
          yhat[,i] <- as.vector(mod[[4]])}
     rownames(monitor) <- Data$mon.names
     ### Warnings
     if(any(is.na(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.")
     if(any(is.nan(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.")
     if(any(is.infinite(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.")
     if(any(!is.finite(Dev)))
          cat("\nWARNING: Deviance has non-finite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat, Deviance=Dev,
          monitor=monitor)
     class(predicted) <- "laplace.ppc"
     return(predicted)
     }

#End
