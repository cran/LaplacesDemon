###########################################################################
# Combine                                                                 #
#                                                                         #
# The purpose of the Combine function is to combine multiple objects of   #
# class demonoid.                                                         #
###########################################################################

Combine <- function(x, Data)
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     len.x <- Thinning <- length(x)
     if(!all(sapply(x, class) == "demonoid")) {
          stop("At least one item in list x is not of class demonoid.")}
     #Note: identical() causes an error below...why???
     #if(!identical(sapply(x, with, Algorithm))) {
     #     cat("\nWarning: Combining different algorithms.\n")}
     Acceptance.Rate <- round(sum(sapply(x, with, Acceptance.Rate) *
               sapply(x, with, Iterations)) /
               sum(sapply(x, with, Iterations)),5)
     Adaptive <- round(mean(sapply(x, with, Adaptive)))
     Algorithm <- x[[1]]$Algorithm
     Covar <- matrix(rowSums(sapply(x, with, Covar) *
               sapply(x, with, Iterations)) /
               sum(sapply(x, with, Iterations)),
               nrow(x[[1]]$Covar), ncol(x[[1]]$Covar))
     Iterations <- sum(sapply(x, with, Iterations))
     Model <- x[[1]]$Model
     Minutes <- sum(sapply(x, with, Minutes))
     Periodicity <- round(mean(sapply(x, with, Periodicity)))
     Status <- round(mean(sapply(x, with, Status)))
     LIV <- x[[1]]$Parameters
     ### Combine
     post <- x[[1]]$Posterior1
     Dev <- x[[1]]$Deviance
     Mon <- x[[1]]$Monitor
     for (i in 2:len.x) {
          post <- rbind(post, x[[i]]$Posterior1)
          Dev <- rbind(Dev, x[[i]]$Deviance)
          Mon <- rbind(Mon, x[[i]]$Monitor)
          }
     ### Thinning
     thinned <- post[seq(from=1, to=nrow(post), by=len.x),]
     Dev <- Dev[seq(from=1, to=length(Dev), by=len.x)]
     Mon <- Mon[seq(from=1, to=nrow(Mon), by=len.x),]
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n")
     burn.start <- trunc(seq(from=1, to=nrow(thinned), by=nrow(thinned)/10))
     geweke <- matrix(9, length(burn.start), LIV)
     geweke.ct <- rep(0, LIV)
     options(warn=-1)
     for (i in 1:length(burn.start))
          {
          thinned2 <- thinned[burn.start[i]:nrow(thinned),]
          test <- try(as.vector(Geweke.Diagnostic(thinned2)), silent=TRUE)
          if(class(test) != "try-error") geweke[i,] <- as.vector(test)
          }
     options(warn=0)
     rm(thinned2)
     geweke <- ifelse(!is.finite(geweke), 9, geweke)
     geweke <- ifelse({geweke > -2} & {geweke < 2}, TRUE, FALSE)
     for (j in 1:LIV) {geweke.ct[j] <- which(geweke[,j] == TRUE)[1]}
     geweke.ct <- ifelse(is.na(geweke.ct), nrow(thinned), geweke.ct)
     BurnIn <- burn.start[max(geweke.ct)]
     BurnIn <- ifelse(is.na(BurnIn), nrow(thinned), BurnIn)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n")
     acf.temp <- matrix(1, trunc(10*log10(nrow(thinned))), LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV)
          {
          temp0 <- acf(thinned[,j], lag.max=NROW(acf.temp), plot=FALSE)
          acf.temp[,j] <- abs(temp0$acf[2:{NROW(acf.temp)+1},,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning
          }
     Rec.Thin <- ifelse(is.na(Rec.Thin), NROW(acf.temp), Rec.Thin)
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     ### Assess ESS for stationary samples
     if(BurnIn < nrow(thinned))
          {
          ESS4 <- ESS(thinned[BurnIn:nrow(thinned),])
          ESS5 <- ESS(Dev[BurnIn:length(Dev)])
          ESS6 <- ESS(Mon[BurnIn:nrow(Mon),])
          }
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Mon <- as.matrix(Mon)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- Summ1[,2] / sqrt(ESS1)
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025))
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500))
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975))
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(Dev)
     Deviance[3] <- sd(Dev) / sqrt(ESS2)
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(Mon[,j])
          Monitor[3] <- sd(Mon[,j]) / sqrt(ESS3[j])
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.5,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data$mon.names[j]
          }
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(BurnIn < nrow(thinned))
          {
          Summ2[,1] <- colMeans(thinned[BurnIn:nrow(thinned),])
          Summ2[,2] <- apply(thinned[BurnIn:nrow(thinned),], 2, sd)
          Summ2[,3] <- Summ2[,2] / sqrt(ESS4)
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned[BurnIn:nrow(thinned),], 2, quantile,
               c(0.025))
          Summ2[,6] <- apply(thinned[BurnIn:nrow(thinned),], 2, quantile,
               c(0.500))
          Summ2[,7] <- apply(thinned[BurnIn:nrow(thinned),], 2, quantile,
               c(0.975))
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev[BurnIn:length(Dev)])
          Deviance[2] <- sd(Dev[BurnIn:length(Dev)])
          Deviance[3] <- sd(Dev[BurnIn:length(Dev)]) / sqrt(ESS5)
          Deviance[4] <- ESS5
          Deviance[5] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[BurnIn:nrow(thinned),j])
               Monitor[2] <- sd(Mon[BurnIn:nrow(thinned),j])
               Monitor[3] <- sd(Mon[BurnIn:nrow(thinned),j]) /
                    sqrt(ESS6[j])
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon[BurnIn:nrow(Mon),j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[BurnIn:nrow(Mon),j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[BurnIn:nrow(Mon),j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]
               }
          }
     ### Column names to samples
     if(ncol(Mon) == length(Data$mon.names)) colnames(Mon) <- Data$mon.names
     if(ncol(thinned) == length(Data$parm.names)) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     if(Algorithm == "Random-Walk Metropolis") {
          cat("Estimating Log of the Marginal Likelihood\n")
          if(BurnIn >= nrow(thinned)) {LML <- LML(Model, Data,
               thinned[nrow(thinned),])}
          else {LML <- LML(Model, Data, Summ2[1:LIV,6])}}
     else {LML <- list(LML=NA, VarCov=NA)}
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=Acceptance.Rate,
          Adaptive=Adaptive,
          Algorithm=Algorithm,
          Call=match.call(),
          Covar=Covar,
          CovarDHis=x[[len.x]]$CovarDHis,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(BurnIn < nrow(thinned)) {
               c(mean(as.vector(Dev[BurnIn:length(Dev)])),
               var(as.vector(Dev[BurnIn:length(Dev)]))/2,
               mean(as.vector(Dev[BurnIn:length(Dev)])) + 
               var(as.vector(Dev[BurnIn:length(Dev)]))/2)}
               else rep(NA,3),
          DR=x[[1]]$DR,
          Initial.Values=x[[len.x]]$Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=Minutes,
          Model=Model,
          Monitor=Mon,
          Parameters=LIV,
          Periodicity=Periodicity,
          Posterior1=thinned,
          Posterior2=thinned[BurnIn:nrow(thinned),],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=nrow(thinned), Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n")
     return(LaplacesDemon.out)
     }

#End
