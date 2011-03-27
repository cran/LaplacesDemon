###########################################################################
# predict.demonoid                                                        #
#                                                                         #
# The purpose of this function is to predict y[new] or y[rep], and later  #
# provide posterior predictive checks for objects of class demonoid.      #
###########################################################################

predict.demonoid <- function(object, Model, Data, ...)
     {
     ### Initial Checks
     if(is.null(object)) stop("The object argument is NULL.\n")
     if(is.null(Model)) stop("The Model argument is NULL.\n")
     if(is.null(Data)) stop("The Data argument is NULL.\n")
     if(is.null(Data$y) & is.null(Data$Y)) stop("Data must have y or Y.\n")
     if(!is.null(Data$y)) y <- as.vector(Data$y)
     if(!is.null(Data$Y)) y <- as.vector(Data$Y)
     ### p(y[rep] | y)
     post <- as.matrix(object$Posterior1)
     if(is.matrix(object$Posterior2) == TRUE) {
          post <- as.matrix(object$Posterior2)}
     yhat <- matrix(NA, length(y), NROW(post))
     for (i in 1:NROW(post)) {
          yhat[,i] <- as.vector(Model(post[i,], Data)[[4]])}
     ### Warnings
     if(sum(is.na(yhat)) > 0) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.")
     if(sum(is.nan(yhat)) > 0) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.")
     if(sum(is.infinite(yhat)) > 0) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat)
     class(predicted) <- "demonoid.ppc"
     return(predicted)
     }

#End
