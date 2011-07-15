###########################################################################
# summary.demonoid.ppc                                                    #
#                                                                         #
# The purpose of this function is to summarize an object of class         #
# demonoid.ppc (posterior predictive check).                              #
###########################################################################

summary.demonoid.ppc <- function(object=NULL, Categorical=FALSE, Rows=NULL,
     Discrep=NULL, d=0, ...)
     {
     if(is.null(object)) stop("The object argument is NULL.\n")
     y <- object$y
     yhat <- object$yhat
     if(is.null(Rows)) Rows <- 1:NROW(y)
     ### Create Continuous Summary Table
     if(Categorical == FALSE) {
          Summ <- matrix(NA, length(y), 8, dimnames=list(1:length(y),
               c("y","Mean","SD","LB","Median","UB","PQ","Discrep")))
          Summ[,1] <- y
          Summ[,2] <- rowMeans(yhat)
          Summ[,3] <- apply(yhat, 1, sd)
          for (i in 1:length(y))
               {
               Summ[i,4] <- quantile(yhat[i,], probs=0.025)
               Summ[i,5] <- quantile(yhat[i,], probs=0.500)
               Summ[i,6] <- quantile(yhat[i,], probs=0.975)
               Summ[i,7] <- mean(yhat[i,] >= y[i])
               }
          ### Discrepancy Statistics
          Concordance <- 1 - mean({{Summ[,7] < 0.025} | {Summ[,7] > 0.975}},
               na.rm=TRUE)
          if(identical(yhat,y)) Concordance <- 1
          Discrepancy.Statistic <- 0
          if(!is.null(Discrep) && {Discrep == "max(yhat[i,]) > max(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- max(yhat[i,]) > max(y)}
               Discrepancy.Statistic <- mean(Summ[,8])}
          if(!is.null(Discrep) && {Discrep == "min(yhat[i,]) < min(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- min(yhat[i,]) < min(y)}
               Discrepancy.Statistic <- mean(Summ[,8])}
          if(!is.null(Discrep) && {Discrep == "round(yhat[i,]) = d"}) {
               for (i in 1:length(y)) {Summ[i,8] <- round(yhat[i,]) == d}
               Discrepancy.Statistic <- mean(Summ[,8])}
          ### Create Output
          Summ.out <- list(Concordance=Concordance,
               Discrepancy.Statistic=round(Discrepancy.Statistic,5),
               Summary=Summ[Rows,])
          cat("Concordance: ", Concordance, "\n")
          cat("Discrepancy Statistic: ", round(Discrepancy.Statistic,5), "\n")
          cat("Records: \n")
          print(Summ[Rows,])
          }
     ### Create Categorical Summary Table
     else {
          catcounts <- table(y)
          sumnames <- rep(NA, length(catcounts)+3)
          sumnames[1] <- "y"
          for (i in 1:length(catcounts)) {
               sumnames[i+1] <- paste("p(yhat=",names(catcounts)[i],")",sep="")}
          sumnames[length(sumnames)-1] <- "Lift"
          sumnames[length(sumnames)] <- "Discrep"
          Summ <- matrix(NA, length(y), length(sumnames),
               dimnames=list(1:length(y), sumnames))
          Summ[,1] <- y
          for (i in 1:length(catcounts)) {
               Summ[,i+1] <- rowSums(yhat == as.numeric(names(catcounts)[i])) /
                    ncol(yhat)}
          Summ[,{ncol(Summ)-1}] <- 1
          for (i in 1:length(y)) {
               Summ[i,{ncol(Summ)-1}] <- Summ[i,
                    grep(Summ[i,1],names(catcounts))+1] / 
                    {as.vector(catcounts[grep(Summ[i,1],names(catcounts))]) /
                    sum(catcounts)} - 1}
          ### Discrepancy Statistics
          Mean.Lift <- mean(Summ[,{ncol(Summ)-1}])
          Discrepancy.Statistic <- 0
          if(!is.null(Discrep) && {Discrep == "p(yhat[i,] != y[i])"}) {
               for (i in 1:length(y)) {Summ[i,ncol(Summ)] <- 1 - 
                    Summ[i, grep(Summ[i,1],names(catcounts))+1]}
               Discrepancy.Statistic <- mean(Summ[,ncol(Summ)])}
          ### Create Output
          Summ.out <- list(Mean.Lift=Mean.Lift,
               Discrepancy.Statistic=round(Discrepancy.Statistic,5),
               Summary=Summ[Rows,])
          cat("Mean Lift: ", Mean.Lift, "\n")
          cat("Discrepancy Statistic: ", round(Discrepancy.Statistic,5), "\n")
          cat("Records: \n")
          print(Summ[Rows,])
          }
     return(invisible(Summ.out))
     }

#End
