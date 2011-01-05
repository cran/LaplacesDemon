###########################################################################
# predict.demonoid                                                        #
#                                                                         #
# The purpose of this function is to predict y[new] or y[rep], and        #
# provide posterior predictive checks.                                    #
###########################################################################

predict.demonoid <- function(object, Model, Data, ...)
     {
     ### p(y[rep] | y)
     post <- as.matrix(object$Posterior2)
     yhat <- matrix(NA,length(Data$y),NROW(post))
     for(i in 1:NROW(post))
          {yhat[,i] <- Model(post[i,], Data)[[4]]}
     ### Create Output
     predicted <- list(y=Data$y, yhat=yhat)
     class(predicted) <- "demonoid.ppc"
     return(predicted)
     }

#End
