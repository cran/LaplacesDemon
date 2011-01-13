###########################################################################
# summary.demonoid.ppc                                                    #
#                                                                         #
# The purpose of this function is to summarize an object of class         #
# demonoid.ppc (posterior predictive check).                              #
###########################################################################

summary.demonoid.ppc <- function(object=NULL, Rows=NULL, ...)
     {
     if(is.null(object)) {cat("ERROR: The object argument is empty.\n")}
     y <- object$y
     yhat <- object$yhat
     if(is.null(Rows)) {Rows <- 1:NROW(y)}
     ### Create Summary Table
     Summ <- matrix(NA, length(y), 7, dimnames=list(1:length(y),
          c("y","Mean","SD","LB","Median","UB","PQ")))
     Summ[,1] <- y
     Summ[,2] <- apply(yhat, 1, mean)
     Summ[,3] <- apply(yhat, 1, sd)
     for(i in 1:length(y))
         {
         Summ[i,4] <- quantile(yhat[i,], probs=0.025)
         Summ[i,5] <- quantile(yhat[i,], probs=0.500)
         Summ[i,6] <- quantile(yhat[i,], probs=0.975)
         Summ[i,7] <- mean(yhat[i,] >= y[i])
       }
     Concordance <- 1 - mean((Summ[,7] < 0.025) | (Summ[,7] > 0.975))
     ### Create Output
     Summ.out <- list(Concordance=Concordance,
          Summary=Summ[Rows,])
     cat("Concordance: ", Concordance, "\n")
     cat("Records: \n")
     print(Summ[Rows,])
     invisible(Summ.out)
     }

#End
