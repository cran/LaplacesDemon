###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of this function is to perform one of four MCMC algorithms  #
# on the log of a joint posterior distribution of a Bayesian model.       #
###########################################################################

LaplacesDemon <- function(Log.Posterior=NULL, Data=NULL, Adaptive=0,
     Covar=NULL, DR=1, Initial.Values=NULL, Iterations=1000,
     Periodicity=1, Status=100, Thinning=1)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     ### Initial Checks
     cat("\nPerforming initial checks...\n")
     if(is.null(Log.Posterior)) {cat("ERROR: A function must be ",
          "entered for Log.Posterior.\n", sep="")}
     if(is.null(Data)) {cat("ERROR: A list containing data must ",
          "be entered for Data.\n", sep="")}
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="")}
     if((Adaptive <= 1) | (Adaptive > Iterations)) {
          Adaptive <- Iterations + 1
          cat("Adaptation will not occur due to the Adaptive argument.\n")}
     if((DR != 0) & (DR != 1)) {
          DR <- 0
          cat("Delayed Rejection will not occur due to the DR argument.\n")}
     if((Periodicity <= 1) | (Periodicity > Iterations)) {
          Periodicity <- Iterations + 1
          cat("Adaptation will not occur due to the Periodicity argument.\n")}
     if((Status < 1) | (Status > Iterations)) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="")}
     if((Thinning < 1) | (Thinning > Iterations)) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="")}
     ### Initial Settings
     Acceptance <- 0
     LP0 <- Log.Posterior(Initial.Values, Data)
     if(is.na(LP0[[1]])) LP0[[1]] <- -999999999
     if(is.infinite(LP0[[2]])) LP0[[2]] <- 999999999
     Dev <- LP0[[2]]
     Mon <- LP0[[3]]
     LIV <- length(Initial.Values)
     post <- matrix(0, Iterations, LIV)
     thinned <- Initial.Values
     post[1,] <- prop <- Initial.Values
     ScaleF <- 2.381204 / sqrt(LIV) #2.381204^2 / LIV
     if(is.matrix(Covar)) {tuning <- diag(Covar); VarCov <- Covar}
     else {
     tuning <- rep(ScaleF, LIV)
     VarCov <- matrix(0, LIV, LIV)
     diag(VarCov) <- ScaleF / LIV}
     Iden.Mat <- VarCov
     ### Determine Algorithm
     if((Adaptive < Iterations) & (DR == 0)) {
          Algorithm <- "Adaptive Metropolis"
          cat("Algorithm: Adaptive Metropolis\n")}
     if((Adaptive < Iterations) & (DR == 1)) {
          Algorithm <- "Delayed Rejection Adaptive Metropolis"
          cat("Algorithm: Delayed Rejection Adaptive Metropolis\n")}
     if((Adaptive >= Iterations) & (DR == 1)) {
          Algorithm <- "Delayed Rejection Metropolis"
          cat("Algorithm: Delayed Rejection Metropolis\n")}
     if((Adaptive >= Iterations) & (DR == 0)) {
          Algorithm <- "Random-Walk Metropolis"
          cat("Algorithm: Random-Walk Metropolis\n")}
     cat("\nLaplace's Demon is beginning to update...\n")
     for (iter in 1:(Iterations-1))
          {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter, sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0)
               {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, LP0[[2]])
               Mon <- rbind(Mon, LP0[[3]])
               }
          ### Prepare next iteration in case of rejection
          post[iter+1,] <- post[iter,]
          ### Log-Posterior of the current state
          LP0 <- Log.Posterior(post[iter,], Data)
          if(is.na(LP0[[1]])) LP0[[1]] <- -999999999
          if(is.infinite(LP0[[2]])) LP0[[2]] <- 999999999
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[1])) {
               if(iter %% Status == 0) {
                   cat(",   Proposal: Multivariate\n")}
               MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov)
               prop <- t(post[iter,] + t(MVNz))}
          if(is.character(MVN.test[1])) {
               if(iter %% Status == 0) {
                    cat(",   Proposal: Independent\n")}
               for (j in 1:LIV) {
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}}
          ### Log-Posterior of the proposed state
          LP1 <- Log.Posterior(prop, Data)
          if(is.na(LP1[[1]])) LP1[[1]] <- LP0[[1]]
          if(is.infinite(LP1[[2]])) LP1[[2]] <- LP0[[1]]
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- LP1[[1]] - LP0[[1]]
          log.alpha <- ifelse(is.na(log.alpha), 0, log.alpha)
          if(log.u < log.alpha)
               {
               post[iter+1,] <- prop
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0)
                    {
                    Dev[length(Dev)] <- LP1[[2]]
                    if(is.vector(Mon)) {Mon[length(Mon)] <- LP1[[3]][1]}
                    if(is.matrix(Mon)) {Mon[NROW(Mon),] <- LP1[[3]]}
                    }
               }
          ### Delayed Rejection: Second Stage Proposals
          if((DR == 1) & (log.u >= log.alpha))
               {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov), silent=TRUE)
               if(is.numeric(MVN.test[1])) {
                    MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov)
                    prop <- t(post[iter,] + t(MVNz))}
               if(is.character(MVN.test[1])) {for (j in 1:LIV) {
                         prop[j] <- rnorm(1, prop[j], tuning[j])}}
               ### Log-Posterior of the proposed state
               LP12 <- Log.Posterior(prop, Data)
               if(is.na(LP12[[1]])) LP12[[1]] <- LP1[[1]]
               if(is.infinite(LP12[[2]])) LP12[[2]] <- LP1[[1]]
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(LP1[[1]] - LP12[[1]]))
               options(warn=0)
               if(is.nan(log.alpha.comp)) log.alpha.comp <- 0
               if(is.infinite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- LP12[[1]] + log.alpha.comp  -
                    (LP0[[1]] + log(1 - exp(LP1[[1]] - LP0[[1]])))
               log.alpha <- ifelse(is.na(log.alpha), 0, log.alpha)
               if(log.u < log.alpha)
                    {
                    post[iter+1,] <- prop
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0)
                         {
                         Dev[length(Dev)] <- LP1[[2]]
                         if(is.vector(Mon)) {
                              Mon[length(Mon)] <- LP1[[3]][1]}
                         if(is.matrix(Mon)) {Mon[NROW(Mon),] <- LP1[[3]]}
                         }
                    }
               }
          ### Shrinkage of Adaptive Proposal Variance
          if((Adaptive < Iterations) & (Acceptance / iter < 0.05))
               {
               VarCov <- VarCov * 0.99
               tuning <- tuning * 0.99
               }
          ### Adapt the Proposal Variance
          if((iter >= Adaptive) & (iter %% Periodicity == 0))
               {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- (ScaleF * cov(post[1:iter,])) +
                    (ScaleF * 1.0E-5 * Iden.Mat)
               ### Univariate Standard Deviations
               for (j in 1:LIV)
                    {
                    tuning[j] <- sqrt(ScaleF * (var(post[1:iter,j])) +
                         ScaleF * 1.0E-5)
                    }
               }
          }
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n")
     burn.start <- trunc(seq(from=1, to=NROW(thinned), by=NROW(thinned)/10))
     geweke <- matrix(9, length(burn.start), LIV)
     geweke.ct <- rep(0, LIV)
     options(warn=-1)
     for (i in 1:length(burn.start))
          {
          thinned2 <- thinned[burn.start[i]:NROW(thinned),]
          test <- try(as.vector(Geweke.Diagnostic(thinned2)), silent=TRUE)
          if(is.numeric(test)) {
               geweke[i,] <- as.vector(Geweke.Diagnostic(thinned2))
          }}
     options(warn=0)
     rm(thinned2)
     geweke <- ifelse(is.na(geweke), 9, geweke)
     geweke <- ifelse(is.nan(geweke), 9, geweke)
     geweke <- ifelse((geweke > -2) & (geweke < 2), TRUE, FALSE)
     for (j in 1:LIV)
          {geweke.ct[j] <- which(geweke[,j] == TRUE)[1]}
     geweke.ct <- ifelse(is.na(geweke.ct), NROW(thinned), geweke.ct)
     BurnIn <- burn.start[max(geweke.ct)]
     BurnIn <- ifelse(is.na(BurnIn), NROW(thinned), BurnIn)
     ### Assess Thinning and Effective Size for all parameter samples
     cat("Assessing Thinning and Effective Size\n")
     acf.temp <- matrix(1, trunc(10*log10(NROW(thinned))), LIV)
     Eff.Size1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV)
          {
          temp0 <- acf(thinned[,j], lag.max=NROW(acf.temp), plot=FALSE)
          acf.temp[,j] <- abs(temp0$acf[2:(NROW(acf.temp)+1),,1])
          Eff.Size1[j] <- Effective.Size(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning
          }
     Rec.Thin <- ifelse(is.na(Rec.Thin), NROW(acf.temp), Rec.Thin)
     ### Assess Effective Size for all deviance and monitor samples
     Eff.Size2 <- Effective.Size(Dev)
     Eff.Size3 <- Effective.Size(Mon)
     ### Assess Effective Size for stationary samples
     if(BurnIn < NROW(thinned))
          {
          Eff.Size4 <- Effective.Size(thinned[BurnIn:NROW(thinned),])
          Eff.Size5 <- Effective.Size(Dev[BurnIn:NROW(Dev),])
          Eff.Size6 <- Effective.Size(Mon[BurnIn:NROW(Mon),])
          }
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     if(is.vector(Mon)) Num.Mon <- 1
     if(is.matrix(Mon)) Num.Mon <- NCOL(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","Eff.Size","LB","Median","UB")))
     Summ1[,1] <- apply(thinned, 2, mean)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- Summ1[,2] / sqrt(Eff.Size1)
     Summ1[,4] <- Eff.Size1
     for (j in 1:LIV) {
          Summ1[j,5] <- as.numeric(quantile(thinned[,j], probs=0.025,
               na.rm=TRUE))
          Summ1[j,6] <- as.numeric(quantile(thinned[,j], probs=0.5,
               na.rm=TRUE))
          Summ1[j,7] <- as.numeric(quantile(thinned[,j], probs=0.975,
               na.rm=TRUE))
          }
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(Dev)
     Deviance[3] <- sd(Dev) / sqrt(Eff.Size2)
     Deviance[4] <- Eff.Size2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     if(Num.Mon == 1)
          {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon)
          Monitor[2] <- sd(Mon)
          Monitor[3] <- sd(Mon) / sqrt(Eff.Size3)
          Monitor[4] <- Eff.Size3
          Monitor[5] <- as.numeric(quantile(Mon, probs=0.025, na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon, probs=0.500, na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon, probs=0.975, na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          }
     if(Num.Mon > 1)
          {
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(Mon[,j])
               Monitor[3] <- sd(Mon[,j]) / sqrt(Eff.Size3[j])
               Monitor[4] <- Eff.Size3[j]
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.5,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ1 <- rbind(Summ1, Monitor)
               }
          }
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","Eff.Size","LB","Median","UB")))
     if(BurnIn < NROW(thinned))
          {
          Summ2[,1] <- apply(thinned[BurnIn:NROW(thinned),], 2, mean)
          Summ2[,2] <- apply(thinned[BurnIn:NROW(thinned),], 2, sd)
          Summ2[,3] <- Summ2[,2] / sqrt(Eff.Size4)
          Summ2[,4] <- Eff.Size4
          for (j in 1:LIV) {
               Summ2[j,5] <- as.numeric(
                    quantile(thinned[BurnIn:NROW(thinned),j], probs=0.025,
                    na.rm=TRUE))
               Summ2[j,6] <- as.numeric(
                    quantile(thinned[BurnIn:NROW(thinned),j], probs=0.500,
                    na.rm=TRUE))
               Summ2[j,7] <- as.numeric(
                    quantile(thinned[BurnIn:NROW(thinned),j], probs=0.975,
                    na.rm=TRUE))
               }
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev[BurnIn:length(Dev)])
          Deviance[2] <- sd(Dev[BurnIn:length(Dev)])
          Deviance[3] <- sd(Dev[BurnIn:length(Dev)]) / sqrt(Eff.Size5)
          Deviance[4] <- Eff.Size5
          Deviance[5] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev[BurnIn:length(Dev)],
               probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          if(Num.Mon == 1)
               {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[BurnIn:NROW(thinned)])
               Monitor[2] <- sd(Mon[BurnIn:NROW(thinned)])
               Monitor[3] <- sd(Mon[BurnIn:NROW(thinned)]) /
                    sqrt(Eff.Size6)
               Monitor[4] <- Eff.Size6
               Monitor[5] <- as.numeric(
                    quantile(Mon[BurnIn:NROW(thinned)], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(
                    quantile(Mon[BurnIn:NROW(thinned)], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(
                    quantile(Mon[BurnIn:NROW(thinned)], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               }
          if(Num.Mon > 1)
               {
               for (j in 1:Num.Mon) {
                    Monitor <- rep(NA,7)
                    Monitor[1] <- mean(Mon[BurnIn:NROW(thinned),j])
                    Monitor[2] <- sd(Mon[BurnIn:NROW(thinned),j])
                    Monitor[3] <- sd(Mon[BurnIn:NROW(thinned),j]) /
                         sqrt(Eff.Size6[j])
                    Monitor[4] <- Eff.Size6[j]
                    Monitor[5] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                         probs=0.025, na.rm=TRUE))
                    Monitor[6] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                         probs=0.500, na.rm=TRUE))
                    Monitor[7] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                         probs=0.975, na.rm=TRUE))
                    Summ2 <- rbind(Summ2, Monitor)
                    }
               }
          }
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,3),
          Adaptive=Adaptive,
          Algorithm=Algorithm,
          Call=match.call(),
          Covar=VarCov,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(BurnIn < NROW(thinned)) {
               c(mean(as.vector(Dev[BurnIn:length(Dev)])),
               var(as.vector(Dev[BurnIn:length(Dev)]))/2,
               mean(as.vector(Dev[BurnIn:length(Dev)])) + 
               var(as.vector(Dev[BurnIn:length(Dev)]))/2)}
               else rep(NA,3),
          DR=DR,
          Iterations=Iterations,
          Log.Posterior=Log.Posterior,
          Minutes=as.vector(time2[3] - time1[3]) / 60,
          Monitor=Mon,
          Parameters=LIV,
          Periodicity=Periodicity,
          Posterior1=thinned,
          Posterior2=thinned[BurnIn:NROW(thinned),],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=NROW(thinned), Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n")
     return(LaplacesDemon.out)
     }

#End
