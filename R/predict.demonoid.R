###########################################################################
# predict.demonoid                                                        #
#                                                                         #
# The purpose of this function is to predict y[new] or y[rep], and        #
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
     if(sum(is.na(object$Posterior2)) == 0)
          {post <- as.matrix(object$Posterior2)}
     else {post <- as.matrix(object$Posterior1)}
     yhat <- matrix(NA, length(y), NROW(post))
     for (i in 1:NROW(post)) {yhat[,i] <- Model(post[i,], Data)[[4]]}
     ### Create Output
     predicted <- list(y=y, yhat=yhat)
     class(predicted) <- "demonoid.ppc"
     return(predicted)
     }

#End
