###########################################################################
# predict.demonoid                                                        #
#                                                                         #
# The purpose of this function is to predict y[rep] and provide posterior #
# predictive checks.                                                      #
###########################################################################

predict.demonoid <- function(object, Log.Posterior, Data, ...)
     {
     ### p(y[rep] | y)
     post <- as.matrix(object$Posterior2)
     yrep <- matrix(NA,length(Data$y),NROW(post))
     for(i in 1:NROW(post))
          {yrep[,i] <- Log.Posterior(post[i,], Data)[[4]]}
     ### Create Summary Table
     Summ <- matrix(NA, length(Data$y), 6, dimnames=list(1:length(Data$y),
          c("Mean","SD","LB","Median","UB","PPPV")))
     Summ[,1] <- apply(yrep, 1, mean)
     Summ[,2] <- apply(yrep, 1, sd)
     for(i in 1:length(Data$y))
         {
         Summ[i,3] <- quantile(yrep[i,], probs=0.025)
         Summ[i,4] <- quantile(yrep[i,], probs=0.500)
         Summ[i,5] <- quantile(yrep[i,], probs=0.975)
         Summ[i,6] <- mean(yrep[i,] >= Data$y[i])
       }
     Concordance <- 1 - mean((Summ[,6] < 0.025) | (Summ[,6] > 0.975))
     ### Create Output
     predicted <- list(Concordance=Concordance, Summary=Summ, yrep=yrep)
     return(predicted)
     }

#End
