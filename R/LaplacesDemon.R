###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of this function is to perform one of four MCMC algorithms  #
# on the logarithm of the unnormalized joint posterior density of a       #
# Bayesian model.                                                         #
###########################################################################

LaplacesDemon <- function(Model=NULL, Data=NULL, Adaptive=0,
     Covar=NULL, DR=1, Initial.Values=NULL, Iterations=1000,
     Periodicity=1, Status=100, Thinning=1)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(is.null(Model)) stop("A function must be entered for Model.\n")
     if(is.null(Data)) stop("A list containing data must be entered for Data.\n")
     if(is.null(Data$mon.names)) stop("In Data, mon.names is NULL.\n")
     if(is.null(Data$parm.names)) stop("In Data, parm.names is NULL.\n")
     if(is.null(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(length(Initial.Values) != length(Data$parm.names)) {
          cat("WARNING: The length of Initial Values differed from ",
               "Data$parm.names.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(is.na(Initial.Values))) {
          cat("WARNING: Initial Values contain missing values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(is.nan(Initial.Values))) {
          cat("WARNING: Initial Values contain NaN values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(is.infinite(Initial.Values))) {
          cat("WARNING: Initial Values contain infinite values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="")}
     if({Adaptive <= 1} || {Adaptive > Iterations}) {
          Adaptive <- Iterations + 1
          cat("Adaptation will not occur due to the Adaptive argument.\n")}
     if({DR != 0} & {DR != 1}) {
          DR <- 0
          cat("Delayed Rejection will not occur due to the DR argument.\n")}
     if({Periodicity <= 1} || {Periodicity > Iterations}) {
          Periodicity <- Iterations + 1
          cat("Adaptation will not occur due to the Periodicity argument.\n")}
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="")}
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="")}
     Mo0 <- Model(Initial.Values, Data)
     if(length(Mo0[[1]]) > 1) stop("Multiple joint posteriors exist!\n")
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x) 
          out <- c(sprintf('"%s" <- ', fname), capture.output(f)) 
          if(grepl("^[<]", tail(out,1))) out <- head( out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are 'vectorized'.\n")
     }
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are 'vectorized'.\n")
     }
     #########################  Initial Settings  #########################
     Acceptance <- 0
     Mo0 <- Model(Initial.Values, Data)
     if(is.na(Mo0[[1]]) || is.nan(Mo0[[1]]) || is.infinite(Mo0[[1]])) {
          Initial.Values <- runif(length(Initial.Values),-3,3)
          Mo0 <- Model(Initial.Values, Data)
          if(is.na(Mo0[[1]]) || is.nan(Mo0[[1]]) || is.infinite(Mo0[[1]])) {
               Initial.Values <- runif(length(Initial.Values),-3,3)
               Mo0 <- Model(Initial.Values, Data)
               }
          }
     if(is.na(Mo0[[1]])) stop("The posterior is a missing value!\n")
     if(is.infinite(Mo0[[1]])) {
          cat("Initial values are being randomized due to an infinite posterior.\n")
          for (i in 1:1000) {
               Initial.Values <- runif(length(Initial.Values),-2,2)
               Mo0 <- Model(Initial.Values, Data)
               if(!is.infinite(Mo0[[1]])) {
                    cat("Initial values now result in a finite posterior.\n")
                    break
                    }
               }
          }
     if(is.infinite(Mo0[[1]])) stop("The posterior is infinite!\n")
     if(is.nan(Mo0[[1]])) stop("The posterior was not a number!\n")
     if(is.na(Mo0[[2]])) stop("The deviance is a missing value!\n")
     if(is.infinite(Mo0[[2]])) stop("The deviance is infinite!\n")
     if(is.nan(Mo0[[2]])) stop("The deviance is not a number!\n")
     if(any(is.na(Mo0[[3]])))
          stop("Monitored variable(s) have a missing value!\n")
     if(any(is.infinite(Mo0[[3]])))
          stop("Monitored variable(s) have an infinite value!\n")
     if(any(is.nan(Mo0[[3]])))
          stop("Monitored variable(s) include a value that is not a number!\n")
     ######################  Laplace Approximation  #######################
     ### Sample Size of Data
     if(!is.null(Data$n)) if(length(Data$n) == 1) N <- Data$n
     if(!is.null(Data$N)) if(length(Data$N) == 1) N <- Data$N
     if(!is.null(Data$y)) N <- nrow(matrix(Data$y))
     if(!is.null(Data$Y)) N <- nrow(matrix(Data$Y))
     if(is.null(N)) stop("Sample size of Data not found in n, N, y, or Y.")
     if({sum(abs(Initial.Values) == 0) == length(Initial.Values)} &
          {N >= 5*length(Initial.Values)}) {
          cat("\nLaplace Approximation will be used on initial values.\n")
          Fit.LA = LaplaceApproximation(Model, Initial.Values, Data)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) * Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n")
          cat("for Laplace's Demon, and the posterior modes are now the initial\n")
          cat("values for Laplace's Demon.\n\n")
          }
     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- Mo0[[2]]
     Mon <- Mo0[[3]]
     LIV <- length(Initial.Values)
     post <- matrix(0, Iterations, LIV)
     thinned <- Initial.Values
     post[1,] <- prop <- Initial.Values
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(is.matrix(Covar)) {tuning <- diag(Covar); VarCov <- Covar}
     else {
          tuning <- rep(ScaleF, LIV)
          VarCov <- matrix(0, LIV, LIV)
          diag(VarCov) <- ScaleF / LIV}
     Iden.Mat <- VarCov
     DiagCovar <- matrix(diag(VarCov), 1, LIV)
     ### Determine Algorithm
     if({Adaptive < Iterations} & {DR == 0}) {
          Algorithm <- "Adaptive Metropolis"
          cat("Algorithm: Adaptive Metropolis\n")}
     else if({Adaptive < Iterations} & {DR == 1}) {
          Algorithm <- "Delayed Rejection Adaptive Metropolis"
          cat("Algorithm: Delayed Rejection Adaptive Metropolis\n")}
     else if({Adaptive >= Iterations} & {DR == 1}) {
          Algorithm <- "Delayed Rejection Metropolis"
          cat("Algorithm: Delayed Rejection Metropolis\n")}
     else {#if({Adaptive >= Iterations} & {DR == 0}) {
          Algorithm <- "Random-Walk Metropolis"
          cat("Algorithm: Random-Walk Metropolis\n")}
     ############################  Begin MCMC  ############################
     cat("\nLaplace's Demon is beginning to update...\n")
     for (iter in 1:{Iterations-1})
          {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter, sep="")}
          ### Save Thinned Samples
          if(iter %% Thinning == 0)
               {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])
               }
          ### Prepare next iteration in case of rejection
          post[iter+1,] <- post[iter,]
          ### Log-Posterior of the current state
          Mo0 <- Model(post[iter,], Data)
          if(is.na(Mo0[[1]])) stop("The posterior is a missing value!\n")
          if(is.infinite(Mo0[[1]])) stop("The posterior is infinite!\n")
          if(is.nan(Mo0[[1]])) stop("The posterior was not a number!\n")
          if(is.na(Mo0[[2]])) stop("The deviance is a missing value!\n")
          if(is.infinite(Mo0[[2]])) stop("The deviance is infinite!\n")
          if(is.nan(Mo0[[2]])) stop("The deviance is not a number!\n")
          if(any(is.na(Mo0[[3]])))
               stop("Monitored variable(s) have a missing value!\n")
          if(any(is.infinite(Mo0[[3]])))
               stop("Monitored variable(s) have an infinite value!\n")
          if(any(is.nan(Mo0[[3]])))
               stop("Monitored variable(s) include a value that is not a number!\n")
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[1]) & ((Acceptance / Iterations) >= 0.05)) {
               if(iter %% Status == 0) {
                   cat(",   Proposal: Multivariate\n")}
               MVNz <- as.vector(MVN.test)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) {
                    cat(",   Proposal: Single-Component\n")}
               prop <- post[iter,]
               j <- round(runif(1,0.5,{LIV+0.49}))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(is.na(Mo1[[1]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(is.infinite(Mo1[[1]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(is.nan(Mo1[[1]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(is.na(Mo1[[2]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(is.infinite(Mo1[[2]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(is.nan(Mo1[[2]])) {Mo1 <- Mo0; prop <- post[iter,]}
          if(any(is.na(Mo1[[3]]))) {Mo1 <- Mo0; prop <- post[iter,]}
          if(any(is.infinite(Mo1[[3]]))) {Mo1 <- Mo0; prop <- post[iter,]}
          if(any(is.nan(Mo1[[3]]))) {Mo1 <- Mo0; prop <- post[iter,]}
          prop <- Mo1[[5]]
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          log.alpha <- ifelse(is.na(log.alpha), 0, log.alpha)
          if(log.u < log.alpha)
               {
               post[iter+1,] <- prop
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0)
                    {
                    Dev[length(Dev)] <- Mo1[[2]]
                    Mon <- as.matrix(Mon)
                    Mon[NROW(Mon),] <- Mo1[[3]]
                    }
               }
          ### Delayed Rejection: Second Stage Proposals
          else if({DR == 1} & {log.u >= log.alpha})
               {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(is.numeric(MVN.test[1]) & ((Acceptance / Iterations) >= 0.05)) {
                    MVNz <- as.vector(MVN.test)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- round(runif(1,0.5,{LIV+0.49}))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo12 <- Model(prop, Data)
               if(is.na(Mo12[[1]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(is.infinite(Mo12[[1]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(is.nan(Mo12[[1]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(is.na(Mo12[[2]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(is.infinite(Mo12[[2]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(is.nan(Mo12[[2]])) {Mo12 <- Mo0; prop <- post[iter,]}
               if(any(is.na(Mo12[[3]]))) {Mo12 <- Mo0; prop <- post[iter,]}
               if(any(is.infinite(Mo12[[3]]))) {Mo12 <- Mo0; prop <- post[iter,]}
               if(any(is.nan(Mo12[[3]]))) {Mo12 <- Mo0; prop <- post[iter,]}
               prop <- Mo12[[5]]
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[[1]] - Mo12[[1]]))
               options(warn=0)
               if(is.nan(log.alpha.comp)) log.alpha.comp <- 0
               if(is.infinite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo12[[1]] + log.alpha.comp  -
                    {Mo0[[1]] + log(1 - exp(Mo1[[1]] - Mo0[[1]]))}
               log.alpha <- ifelse(is.na(log.alpha), 0, log.alpha)
               if(log.u < log.alpha)
                    {
                    post[iter+1,] <- prop
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0)
                         {
                         Dev[length(Dev)] <- Mo1[[2]]
                         Mon <- as.matrix(Mon)
                         Mon[NROW(Mon),] <- Mo1[[3]]
                         }
                    }
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05})
               {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}
               }
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0})
               {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV)
                    {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)
                    }
               }
          }
     #########################  MCMC is Finished  #########################
     ### Warnings (After Updating)
     if(Acceptance == 0) cat("\nWARNING: All proposals were rejected.\n")
     ### Real Values
     thinned <- ifelse(is.na(thinned), 0, thinned)
     thinned <- ifelse(is.infinite(thinned), 0, thinned)
     thinned <- ifelse(is.nan(thinned), 0, thinned)
     Dev <- ifelse(is.na(Dev), 0, Dev)
     Dev <- ifelse(is.infinite(Dev), 0, Dev)
     Dev <- ifelse(is.nan(Dev), 0, Dev)
     Mon <- ifelse(is.na(Mon), 0, Mon)
     Mon <- ifelse(is.infinite(Mon), 0, Mon)
     Mon <- ifelse(is.nan(Mon), 0, Mon)
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
          if(is.numeric(test)) geweke[i,] <- as.vector(test)
          }
     options(warn=0)
     rm(thinned2)
     geweke <- ifelse(is.na(geweke), 9, geweke)
     geweke <- ifelse(is.nan(geweke), 9, geweke)
     geweke <- ifelse({geweke > -2} & {geweke < 2}, TRUE, FALSE)
     for (j in 1:LIV) {geweke.ct[j] <- which(geweke[,j] == TRUE)[1]}
     geweke.ct <- ifelse(is.na(geweke.ct), NROW(thinned), geweke.ct)
     BurnIn <- burn.start[max(geweke.ct)]
     BurnIn <- ifelse(is.na(BurnIn), NROW(thinned), BurnIn)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n")
     acf.temp <- matrix(1, trunc(10*log10(NROW(thinned))), LIV)
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
     if(BurnIn < NROW(thinned))
          {
          ESS4 <- ESS(thinned[BurnIn:NROW(thinned),])
          ESS5 <- ESS(Dev[BurnIn:NROW(Dev),])
          ESS6 <- ESS(Mon[BurnIn:NROW(Mon),])
          }
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Mon <- as.matrix(Mon)
     Num.Mon <- NCOL(Mon)
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
     if(BurnIn < NROW(thinned))
          {
          Summ2[,1] <- colMeans(thinned[BurnIn:NROW(thinned),])
          Summ2[,2] <- apply(thinned[BurnIn:NROW(thinned),], 2, sd)
          Summ2[,3] <- Summ2[,2] / sqrt(ESS4)
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned[BurnIn:NROW(thinned),], 2, quantile,
               c(0.025))
          Summ2[,6] <- apply(thinned[BurnIn:NROW(thinned),], 2, quantile,
               c(0.500))
          Summ2[,7] <- apply(thinned[BurnIn:NROW(thinned),], 2, quantile,
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
               Monitor[1] <- mean(Mon[BurnIn:NROW(thinned),j])
               Monitor[2] <- sd(Mon[BurnIn:NROW(thinned),j])
               Monitor[3] <- sd(Mon[BurnIn:NROW(thinned),j]) /
                    sqrt(ESS6[j])
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[BurnIn:NROW(Mon),j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]
               }
          }
     ### Column names to samples
     if(NCOL(Mon) == length(Data$mon.names)) colnames(Mon) <- Data$mon.names
     if(NCOL(thinned) == length(Data$parm.names)) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     if(Algorithm == "Random-Walk Metropolis") {
          cat("Estimating Log of the Marginal Likelihood\n")
          if(BurnIn >= NROW(thinned)) {LML <- LML(Model, Data,
               thinned[nrow(thinned),])}
          else {LML <- LML(Model, Data, Summ2[1:LIV,6])}}
     else {LML <- list(LML=NA, VarCov=NA)}
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,7),
          Adaptive=Adaptive,
          Algorithm=Algorithm,
          Call=match.call(),
          Covar=VarCov,
          CovarDHis=DiagCovar,
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
          Initial.Values=Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Model=Model,
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
