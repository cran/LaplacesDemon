###########################################################################
# predict.pmc                                                             #
#                                                                         #
# The purpose of the predict.pmc function is to predict y[new] or y[rep], #
# and later provide posterior predictive checks for objects of class pmc. #
###########################################################################

predict.pmc <- function(object, Model, Data, ...)
     {
     ### Initial Checks
     if(missing(object)) stop("The object argument is required.")
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data$y) & is.null(Data$Y)) stop("Data must have y or Y.")
     if(!is.null(Data$y)) y <- as.vector(Data$y)
     if(!is.null(Data$Y)) y <- as.vector(Data$Y)
     ### p(y[rep] | y)
     post <- as.matrix(object$Posterior2)
     Dev <- rep(NA, nrow(post))
     yhat <- matrix(NA, length(y), nrow(post))
     lengthcomp <- as.vector(Model(post[1,], Data)[[4]])
     if(!identical(length(lengthcomp), length(y)))
          stop("y and yhat differ in length.")
     for (i in 1:nrow(post)) {
          mod <- Model(post[i,], Data)
          Dev[i] <- as.vector(mod[[2]])
          yhat[,i] <- as.vector(mod[[4]])}
     ### Warnings
     if(any(is.na(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.na(yhat)), " missing values.\n")
     if(any(is.nan(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.nan(yhat)), " non-numeric (NaN) values.\n")
     if(any(is.infinite(yhat))) cat("\nWARNING: Output matrix yhat has ",
          sum(is.infinite(yhat)), " infinite values.\n")
     if(any(!is.finite(Dev)))
          cat("\nWARNING: Deviance has non-finite values.")
     ### Create Output
     predicted <- list(y=y, yhat=yhat, Deviance=Dev)
     class(predicted) <- "pmc.ppc"
     return(predicted)
     }

#End
