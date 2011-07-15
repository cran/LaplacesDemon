###########################################################################
# predict.laplace                                                         #
#                                                                         #
# The purpose of this function is to predict y[new] or y[rep], and later  #
# provide posterior predictive checks for objects of class laplace.       #
###########################################################################

predict.laplace <- function(object, Model, Data, Samples=1000, ...)
     {
     ### Initial Checks
     if(is.null(object)) stop("The object argument is NULL.\n")
     if(is.null(Model)) stop("The Model argument is NULL.\n")
     if(is.null(Data)) stop("The Data argument is NULL.\n")
     if(is.null(Data$y) & is.null(Data$Y)) stop("Data must have y or Y.\n")
     if(!is.null(Data$y)) y <- as.vector(Data$y)
     if(!is.null(Data$Y)) y <- as.vector(Data$Y)
     ### Create Posterior Samples of Parameters
     post <- matrix(0, Samples, length(object$Initial.Values))
     for (i in 1:length(object$Initial.Values)) {
          post[,i] <- rnorm(Samples, object$Summary[i,1],
               sqrt(diag(object$Covar)[i]))
          post[,i] <- ifelse(is.na(post[,i]), object$Summary[i,1],
               post[,i])
          post[,i] <- ifelse(is.nan(post[,i]), object$Summary[i,1],
               post[,i])
          post[,i] <- ifelse(is.infinite(post[,i]), object$Summary[i,1],
               post[,i])
          }
     ### p(y[rep] | y), Deviance, and Monitors
     deviance <- rep(NA, Samples)
     monitor <- matrix(NA, length(Data$mon.names), Samples)
     lengthcomp <- as.vector(Model(post[1,], Data)[[4]])
     if(length(lengthcomp) != length(y)) stop("y and yhat differ in length.\n")
     yhat <- matrix(NA, length(y), Samples)
     for (i in 1:Samples) {
          temp <- Model(post[i,], Data)
          deviance[i] <- as.vector(temp[[2]])
          monitor[,i] <- as.vector(temp[[3]])
          yhat[,i] <- as.vector(temp[[4]])
          }
     rownames(monitor) <- Data$mon.names
     ### Warnings
     if(any(is.na(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.")
     if(any(is.nan(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.")
     if(any(is.infinite(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat, deviance=deviance,
          monitor=monitor)
     class(predicted) <- "laplace.ppc"
     return(predicted)
     }

#End
