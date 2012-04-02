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
     if(missing(x)) stop("The x argument is required.")
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
          if(is.matrix(x$Posterior2) == FALSE) {
               post <- cbind(x$Posterior1, x$Monitor)
               cat("\nWARNING: Non-stationary samples used.\n\n")}
          else {
               post <- cbind(x$Posterior2,
                    x$Monitor[x$Rec.BurnIn.Thinned:nrow(x$Monitor),])}
          colnames(post) <- c(colnames(x$Posterior1), colnames(x$Monitor))
          ### Selecting Parms
          if(is.null(Parms)) {keepcols <- 1:ncol(post)}
          else {
               Parms <- sub("\\[","\\\\[",Parms)
               Parms <- sub("\\]","\\\\]",Parms)
               Parms <- sub("\\.","\\\\.",Parms)
               if(length(grep(Parms[1], colnames(post))) == 0)
                    stop("Parameter in Parms does not exist.")
               keepcols <- grep(Parms[1], colnames(post))
               if(length(Parms) > 1) {
                    for (i in 2:length(Parms)) {
                         if(length(grep(Parms[i], colnames(post))) == 0)
                              stop("Parameter in Parms does not exist.")
                         keepcols <- c(keepcols,
                              grep(Parms[i], colnames(post)))}}}
          temp <- colnames(post)[keepcols]
          post <- post[,keepcols]
          colnames(post) <- temp
          #Correlation Table
          options(warn=-1); postcor <- cor(post); options(warn=0)
          #Summary Table
          Summ <- matrix(NA, ncol(post), 6)
          rownames(Summ) <- colnames(post)
          colnames(Summ) <- c("p(theta > 0)","N.Modes","Kurtosis",
               "Skewness","Burn-In","IAT")
          options(warn=-1)
          for (i in 1:ncol(post)) {
               Summ[i,1] <- mean(post[,i] > 0)
               Summ[i,2] <- length(Modes(post[,i])[[1]])
               Summ[i,3] <- round(kurtosis(post[,i]),3)
               Summ[i,4] <- round(skewness(post[,i]),3)
               Summ[i,6] <- round(IAT(post[,i]),3)}
          Summ[,5] <- burnin(post)
          options(warn=0)
          }
     ### class laplace
     if(class(x) == "laplace") {
          ### Posterior
          if(any(is.na(x$Posterior)))
               stop("Posterior samples do not exist.")
          post <- x$Summary1
          ### Selecting Parms
          if(is.null(Parms)) {keeprows <- 1:nrow(post)}
          else {
               Parms <- sub("\\[","\\\\[",Parms)
               Parms <- sub("\\]","\\\\]",Parms)
               Parms <- sub("\\.","\\\\.",Parms)
               if(length(grep(Parms[1], rownames(post))) == 0)
                    stop("Parameter in Parms does not exist.")
               keeprows <- grep(Parms[1], rownames(post))
               if(length(Parms) > 1) {
                    for (i in 2:length(Parms)) {
                         if(length(grep(Parms[i], rownames(post))) == 0)
                              stop("Parameter in Parms does not exist.")
                         keeprows <- c(keeprows,
                              grep(Parms[i], rownames(post)))}}}
          temp <- rownames(post)[keeprows]
          post <- post[keeprows,]
          rownames(post) <- temp
          Posterior <- x$Posterior
          colnames(Posterior) <- rownames(post)
          ### Correlation Table
          options(warn=-1); postcor <- cor(Posterior); options(warn=0)
          ### Summary Table
          Summ <- matrix(NA, nrow(post), 6)
          rownames(Summ) <- rownames(post)
          colnames(Summ) <- c("p(theta > 0)","N.Modes","Kurtosis",
               "Skewness","Burn-In","IAT")
          options(warn=-1)
          for (i in 1:ncol(Posterior)) {
               Summ[i,1] <- mean(Posterior[,i] > 0)
               Summ[i,2] <- length(Modes(Posterior[,i])[[1]])
               Summ[i,3] <- round(kurtosis(Posterior[,i]),3)
               Summ[i,4] <- round(skewness(Posterior[,i]),3)
               Summ[i,5] <- 0
               Summ[i,6] <- round(IAT(Posterior[,i]),3)}
          options(warn=0)
          }
     #Output
     out <- list(Posterior.Correlation=postcor, Posterior.Summary=Summ)
     class(out) <- "posteriorchecks"
     return(out)
     }

#End
