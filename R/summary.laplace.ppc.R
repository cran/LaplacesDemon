###########################################################################
# summary.laplace.ppc                                                     #
#                                                                         #
# The purpose of this function is to summarize an object of class         #
# laplace.ppc (posterior predictive check).                               #
###########################################################################

summary.laplace.ppc <- function(object=NULL, Rows=NULL, Discrep=NULL,
     d=0, ...)
     {
     if(is.null(object)) stop("The object argument is NULL.\n")
     y <- object$y
     yhat <- object$yhat
     deviance <- object$deviance
     monitor <- object$monitor
     if(is.null(Rows)) {Rows <- 1:NROW(y)}
     ### Create Summary Table for y and yhat
     Summ <- matrix(NA, length(y), 8, dimnames=list(1:length(y),
          c("y","Mean","SD","LB","Median","UB","PQ","Discrep")))
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
     ### Discrepancy Statistics
     Concordance <- 1 - mean(((Summ[,7] < 0.025) | (Summ[,7] > 0.975)),
          na.rm=TRUE)
     Discrepancy.Statistic <- 0
     if(!is.null(Discrep) && (Discrep == "max(yhat[i,]) > max(y)")) {
          for (i in 1:length(y)) {Summ[i,8] <- max(yhat[i,]) > max(y)}
          Discrepancy.Statistic <- mean(Summ[,8])}
     if(!is.null(Discrep) && (Discrep == "min(yhat[i,]) < min(y)")) {
          for (i in 1:length(y)) {Summ[i,8] <- min(yhat[i,]) < min(y)}
          Discrepancy.Statistic <- mean(Summ[,8])}
     if(!is.null(Discrep) && (Discrep == "round(yhat[i,]) = d")) {
          for (i in 1:length(y)) {Summ[i,8] <- round(yhat[i,]) == d}
          Discrepancy.Statistic <- mean(Summ[,8])}
     ### Deviance
     dev <- mean(deviance)
     pD <- var(deviance) / 2
     DIC <- dev + pD
     ### Create Summary Table for monitored variables
     Mon <- matrix(NA, nrow(monitor), 5,
          dimnames=list(c(rownames(monitor)),
          c("Mean","SD","LB","Median","UB")))
     for (i in 1:nrow(monitor)) {
          Mon[i,1] <- mean(monitor[i,])
          Mon[i,2] <- sd(monitor[i,])
          Mon[i,3] <- quantile(monitor[i,], probs=0.025)
          Mon[i,4] <- quantile(monitor[i,], probs=0.500)
          Mon[i,5] <- quantile(monitor[i,], probs=0.975)
          }
     ### Create Output
     Summ.out <- list(Concordance=Concordance,
          Discrepancy.Statistic=round(Discrepancy.Statistic,5),
          Summary=Summ[Rows,])
     cat("Concordance: ", Concordance, "\n")
     cat("DIC (Dbar): ", round(dev,3), "\n")
     cat("DIC (pD): ", round(pD,3), "\n")
     cat("DIC (DIC): ", round(DIC,3), "\n")
     cat("Discrepancy Statistic: ", round(Discrepancy.Statistic,5), "\n")
     cat("Monitors:\n")
     print(Mon)
     cat("\n\nRecords:\n")
     print(Summ[Rows,])
     invisible(Summ.out)
     }

#End
