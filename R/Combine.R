###########################################################################
# Combine                                                                 #
#                                                                         #
# The purpose of the Combine function is to combine multiple objects of   #
# class demonoid.                                                         #
###########################################################################

Combine <- function(x, Data, Thinning=1)
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     if(missing(Data)) stop("Data is a required argument.")
     Thinning <- abs(round(Thinning))
     len.x <- length(x)
     if(!all(sapply(x, class) == "demonoid")) {
          stop("At least one item in list x is not of class demonoid.")}
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
     Minutes <- max(sapply(x, with, Minutes))
     Periodicity <- round(mean(sapply(x, with, Periodicity)))
     Status <- round(mean(sapply(x, with, Status)))
     LIV <- x[[1]]$Parameters
     ### Combine
     thinned <- x[[1]]$Posterior1
     Dev <- x[[1]]$Deviance
     Mon <- x[[1]]$Monitor
     for (i in 2:len.x) {
          thinned <- rbind(thinned, x[[i]]$Posterior1)
          Dev <- rbind(Dev, x[[i]]$Deviance)
          Mon <- rbind(Mon, x[[i]]$Monitor)
          }
     ### Thinning
     if(Thinning > 1) {
          thinned <- Thin(thinned, By=Thinning)
          Dev <- Thin(Dev, By=Thinning)
          Mon <- Thin(Mon, By=Thinning)}
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n")
     burn.start <- trunc(seq(from=1, to=nrow(thinned),
          by=nrow(thinned)/10))
     geweke <- matrix(9, length(burn.start), LIV)
     geweke.ct <- rep(0, LIV)
     options(warn=-1)
     for (i in 1:length(burn.start)) {
          thinned2 <- matrix(thinned[burn.start[i]:nrow(thinned),],
               nrow(thinned)-burn.start[i]+1, ncol(thinned))
          test <- try(as.vector(Geweke.Diagnostic(thinned2)), silent=TRUE)
          if(!inherits(test, "try-error")) geweke[i,] <- as.vector(test)}
     options(warn=0)
     rm(thinned2)
     geweke <- ifelse(!is.finite(geweke), 9, geweke)
     geweke <- abs(geweke) < 2
     for (j in 1:LIV) {geweke.ct[j] <- which(geweke[,j] == TRUE)[1]}
     geweke.ct <- ifelse(is.na(geweke.ct), nrow(thinned), geweke.ct)
     BurnIn <- burn.start[max(geweke.ct)]
     BurnIn <- ifelse(is.na(BurnIn), nrow(thinned), BurnIn)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n")
     acf.temp <- matrix(1, trunc(10*log10(nrow(thinned))), LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=nrow(acf.temp), plot=FALSE)
          acf.temp[,j] <- abs(temp0$acf[2:{nrow(acf.temp)+1},,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin <- ifelse(is.na(Rec.Thin), nrow(acf.temp), Rec.Thin)
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     ### Assess ESS for stationary samples
     if(BurnIn < nrow(thinned)) {
          ESS4 <- ESS(thinned[BurnIn:nrow(thinned),])
          ESS5 <- ESS(Dev[BurnIn:nrow(thinned),])
          ESS6 <- ESS(Mon[BurnIn:nrow(thinned),]) }
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Mon <- as.matrix(Mon)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025))
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500))
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975))
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(Dev)
     temp <- try(MCSE(as.vector(Dev)), silent=TRUE)
     if(inherits(temp, "try-error"))
          temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(Mon[,j])
          temp <- try(MCSE(Mon[,j]), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data$mon.names[j]
          }
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(BurnIn < nrow(thinned)) {
          thinned2 <- matrix(thinned[BurnIn:nrow(thinned),],
               nrow(thinned)-BurnIn+1, ncol(thinned))
          Dev2 <- matrix(Dev[BurnIn:nrow(thinned),],
               nrow(thinned)-BurnIn+1, ncol(Dev))
          Mon2 <- matrix(Mon[BurnIn:nrow(thinned),],
               nrow(thinned)-BurnIn+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- apply(thinned2, 2, sd)
          Summ2[,3] <- 0
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i],
                    method="sample.variance")}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(as.vector(Dev2), method="sample.variance")
          Deviance[3] <- temp
          Deviance[4] <- ESS5
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025,
               na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500,
               na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975,
               na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
               if(inherits(temp, "try-error"))
                    temp <- MCSE(as.vector(Mon[,j]),
                         method="sample.variance")
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]}
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data$mon.names)))
          colnames(Mon) <- Data$mon.names
     if(identical(ncol(thinned), length(Data$parm.names))) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(({Algorithm == "Delayed Rejection Metropolis"} |
          {Algorithm == "Metropolis-within-Gibbs"} | 
          {Algorithm == "Random-Walk Metropolis"} |
          {Algorithm == "Sequential Metropolis-within-Gibbs"} |
          {Algorithm == "t-walk"}) &
          {BurnIn < nrow(thinned)}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(theta=thinned2,
               LL=as.vector(Dev2)*(-1/2),
               method="NSIS")}
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance.Rate/Iterations,7),
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
               c(mean(as.vector(Dev2)),
               var(as.vector(Dev2))/2,
               mean(as.vector(Dev2)) + 
               var(as.vector(Dev2))/2)}
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
