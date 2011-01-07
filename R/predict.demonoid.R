###########################################################################
# predict.demonoid                                                        #
#                                                                         #
# The purpose of this function is to predict y[new] or y[rep], and        #
# provide posterior predictive checks.                                    #
###########################################################################

predict.demonoid <- function(object, Model, Data, ...)
     {
     ### Initial Checks
     if(is.null(object)) {cat("ERROR: The object argument is null.\n")}
     if(is.null(Model)) {cat("ERROR: The Model argument is null.\n")}
     if(is.null(Data)) {cat("ERROR: The Data argument is null.\n")}
     ### p(y[rep] | y)
     if(sum(is.na(object$Posterior2)) == 0)
          {post <- as.matrix(object$Posterior2)}
     else {post <- as.matrix(object$Posterior1)}
     yhat <- matrix(NA,length(Data$y),NROW(post))
     for(i in 1:NROW(post))
          {yhat[,i] <- Model(post[i,], Data)[[4]]}
     ### Create Output
     predicted <- list(y=Data$y, yhat=yhat)
     class(predicted) <- "demonoid.ppc"
     return(predicted)
     }

#End
