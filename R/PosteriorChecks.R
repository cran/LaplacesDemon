###########################################################################
# PosteriorChecks                                                         #
#                                                                         #
# The purpose of the PosteriorChecks function is to provide additional    #
# checks of the posterior, including the probability that each theta is   #
# greater than zero, kurtosis, and skewness. This function requires an    #
# object of class demonoid or laplace. The kurtosis and skewness checks   #
# return nothing (NA) for class laplace, because parameters from Laplace  #
# Approximation are normally distributed, by definition.                  #
###########################################################################

PosteriorChecks <- function(x, Parms=NULL)
     {
     ### Initial Checks
     if(is.null(x)) stop("The x argument is NULL.\n")
     if((class(x) != "demonoid") & (class(x) != "laplace"))
          stop("An object of class demonoid or laplace is required.")
     #Kurtosis and Skewness Functions
     kurtosis <- function(x) {  
          m4 <- mean((x-mean(x))^4) 
          kurt <- m4/(sd(x)^4)-3  
          return(kurt)}
     skewness <-  function(x) {
          m3 <- mean((x-mean(x))^3)
          skew <- m3/(sd(x)^3)
          return(skew)}
     ### class demonoid
     if(class(x) == "demonoid") {
          ### Posterior and Monitors
          if(x$Thinned.Samples >= nrow(x$Posterior1)) {
               post <- cbind(x$Posterior1, x$Monitor)
               cat("\nWARNING: Non-stationary samples used.\n\n")}
          if(x$Thinned.Samples < nrow(x$Posterior1)) {
               post <- cbind(x$Posterior1[x$Thinned.Samples:nrow(x$Posterior1),],
                    x$Monitor[x$Thinned.Samples:nrow(x$Posterior1),])}
          colnames(post) <- c(colnames(x$Posterior1), colnames(x$Monitor))
          ### Selecting Parms
          if(is.null(Parms)) {keepcols <- 1:ncol(post)}
          else {
               Parms <- sub("\\[","\\\\[",Parms)
               Parms <- sub("\\]","\\\\]",Parms)
               Parms <- sub("\\.","\\\\.",Parms)
               keepcols <- grep(Parms[1], colnames(post))
               if(length(Parms) > 1) {
                    for (i in 2:length(Parms)) {
                         keepcols <- c(keepcols,
                              grep(Parms[i], colnames(post)))}}}
          temp <- colnames(post)[keepcols]
          post <- post[,keepcols]
          colnames(post) <- temp
          #Correlation Table
          postcor <- cor(post)
          #Summary Table
          Summ <- matrix(NA, ncol(post), 3)
          rownames(Summ) <- colnames(post)
          colnames(Summ) <- c("p(theta > 0)","Kurtosis","Skewness")
          for (i in 1:ncol(post)) {
               Summ[i,1] <- mean(post[,i] > 0)
               Summ[i,2] <- round(kurtosis(post[,i]),3)
               Summ[i,3] <- round(skewness(post[,i]),3)}
          }
     ### class laplace
     if(class(x) == "laplace") {
          ### Posterior
          post <- x$Summary
          ### Selecting Parms
          if(is.null(Parms)) {keeprows <- 1:nrow(post)}
          else {
               Parms <- sub("\\[","\\\\[",Parms)
               Parms <- sub("\\]","\\\\]",Parms)
               Parms <- sub("\\.","\\\\.",Parms)
               keeprows <- grep(Parms[1], rownames(post))
               if(length(Parms) > 1) {
                    for (i in 2:length(Parms)) {
                         keeprows <- c(keeprows,
                              grep(Parms[i], rownames(post)))}}}
          temp <- rownames(post)[keeprows]
          post <- post[keeprows,]
          rownames(post) <- temp
          ### Correlation Table
          postcor <- NA
          ### Summary Table
          Summ <- matrix(NA, nrow(post), 3)
          rownames(Summ) <- rownames(post)
          colnames(Summ) <- c("p(theta > 0)","Kurtosis","Skewness")
          Summ[,1] <- round(pnorm(0, mean=post[,1], sd=post[,2],
               lower.tail=FALSE),3)
          }
     #Output
     out <- list(Posterior.Correlation=postcor, Posterior.Summary=Summ)
     return(out)
     }

#End
