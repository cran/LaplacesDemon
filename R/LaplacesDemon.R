###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=100000, Status=1000, Thinning=100, Algorithm="RWM",
     Specs=NULL)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     LDcall <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(missing(Model)) stop("A function must be entered for Model.")
     if(missing(Data))
          stop("A list containing data must be entered for Data.")
     if(is.null(Data$mon.names)) stop("In Data, mon.names is NULL.")
     if(is.null(Data$parm.names)) stop("In Data, parm.names is NULL.")
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     if(missing(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(!identical(length(Initial.Values), length(Data$parm.names))) {
          cat("WARNING: The length of Initial Values differed from ",
               "Data$parm.names.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     Iterations <- round(Iterations)
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="")}
     Status <- round(Status)
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="")}
     Thinning <- round(Thinning)
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="")}
     if(Algorithm %in% c("AM","AMM","AMWG","DRAM","DRM","MWG","RAM",
          "RWM","SAMWG","SMWG","twalk","USAMWG","USMWG")) {
          if(Algorithm == "AM") {
               Algorithm <- "Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("Adaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("Periodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("Adaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("Periodicity was misspecified and changed to 100.\n")}
               if(Specs[["w"]] == Specs[[3]]) {
                    w <- Specs[[3]]
                    if(w <= 0 || w >= 1) w <- 0.05}
               else {
                    w <- 0.05
                    cat("w was misspecified and changed to 0.05.\n")}}
          if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[1]]) Periodicity <- Specs[[1]]
               else {
                    Periodicity <- 50
                    cat("Periodicity was misspecified and changed to 50.\n")}}
          if(Algorithm == "DRAM") {
               Algorithm <- "Delayed Rejection Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("Adaptive was misspecified and changed to 20.\n")}
               DR <- 1
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("Periodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "DRM") {
               Algorithm <- "Delayed Rejection Metropolis"
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               if(Specs[["alpha.star"]] == Specs[[1]]) {
                    alpha.star <- Specs[[1]]
                    if(alpha.star <= 0 || alpha.star >= 1) {
                         cat("alpha.star not in (0,1). Changed to 0.234.\n")
                         alpha.star <- 0.234}
                    if({length(alpha.star) != 1} &
                         {length(alpha.star) != length(Initial.Values)}) {
                         cat("Length of alpha.star is wrong. Changed to 1.\n")
                         alpha.star <- alpha.star[1]}}
               else {
                    alpha.star <- 0.234
                    cat("alpha.star was misspecified and changed to 0.234.\n")}
               Adaptive <- 2
               DR <- 0
               if(Specs[["Dist"]] == Specs[[2]]) {
                    Dist <- Specs[[2]]
                    if(Dist != "t" & Dist != "N") {
                         cat("Dist was not t or N, and changed to N.\n")
                         Dist <- "N"}}
               else {
                    Dist <- "N"
                    cat("Dist was not t or N, and changed to N.\n")
                    }
               if(Specs[["gamma"]] == Specs[[3]]) {
                    gamma <- Specs[[3]]
                    if(gamma <= 0.5 || gamma > 1) {
                         cat("gamma not in (0.5,1]. Changed to 0.66.\n")
                         gamma <- 0.66}}
               else {
                    gamma <- 0.66
                    cat("gamma was misspecified and changed to 0.66.\n")}
               if(Specs[["Periodicity"]] == Specs[[4]]) Periodicity <- Specs[[4]]
               else {
                    Periodicity <- 10
                    cat("Periodicity was misspecified and changed to 10.\n")}}
          if(Algorithm == "RWM") {
               Algorithm <- "Random-Walk Metropolis"
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "SAMWG") {
               Algorithm <- "Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 50
                    cat("Periodicity was misspecified and changed to 50.\n")}}
          if(Algorithm == "SMWG") {
               Algorithm <- "Sequential Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)}
          if(Algorithm == "twalk") {
               Algorithm <- "t-walk"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(Specs[["n1"]] == Specs[[1]]) {
                    n1 <- Specs[[1]]
                    if(n1 < 1) {
                         cat("n1 must be at least 1. Changed to 4.\n")
                         n1 <- 4}}
               else {
                    n1 <- 4
                    cat("n1 was misspecified and changed to 4.\n")}
               if(Specs[["at"]] == Specs[[2]]) {
                    at <- Specs[[2]]
                    if(at <= 0) {
                         cat("at must be positive. Changed to 6.\n")
                         at <- 6}}
               else {
                    at <- 6
                    cat("at was misspecified and changed to 6.\n")}
               if(Specs[["aw"]] == Specs[[3]]) {
                    aw <- Specs[[3]]
                    if(aw <= 0) {
                         cat("aw must be positive. Changed to 1.5.\n")
                         at <- 1.5}}
               else {
                    aw <- 1.5
                    cat("aw was misspecified and changed to 1.5.\n")}
               Periodicity <- 1}
          if(Algorithm == "USAMWG") {
               Algorithm <- "Updating Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 50
                    cat("Periodicity was misspecified and changed to 50.\n")}
               if({class(Specs[["Fit"]]) == "demonoid"} &
                    {class(Specs[[3]]) == "demonoid"}) Fit <- Specs[[3]]
               if(Specs[["Begin"]] == Specs[[4]]) Begin <- Specs[[4]]}
          if(Algorithm == "USMWG") {
               Algorithm <- "Updating Sequential Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               Periodicity <- Iterations + 1
               if(all(Specs[["Dyn"]] == Specs[[1]])) Dyn <- Specs[[1]]
               else stop("Dyn was misspecified.")
               if(!is.matrix(Dyn)) Dyn <- as.matrix(Dyn)
               if({class(Specs[["Fit"]]) == "demonoid"} &
                    {class(Specs[[2]]) == "demonoid"}) Fit <- Specs[[2]]
               if(Specs[["Begin"]] == Specs[[3]]) Begin <- Specs[[3]]}
               }
     else {cat("Unknown algorithm has been changed to Random-Walk Metropolis.\n")
          Algorithm <- "Random-Walk Metropolis"
          Adaptive <- Iterations + 1
          DR <- 0
          Periodicity <- Iterations + 1
          }
     Adaptive <- round(Adaptive)
     if({Adaptive < 1} || {Adaptive > Iterations}) 
          Adaptive <- Iterations + 1
     Periodicity <- round(Periodicity)
     if({Periodicity < 1} || {Periodicity > Iterations}) 
          Periodicity <- Iterations + 1
     Mo0 <- Model(Initial.Values, Data)
     if(length(Mo0[[1]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(Mo0[[3]]), length(Data$mon.names))) {
          stop("Length of mon.names differs from length of monitors.")}
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out, 1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are 'vectorized'.\n")}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are 'vectorized'.\n")}
     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[[1]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          Initial.Values <- GIV(Model, Data)
          Mo0 <- Model(Initial.Values, Data)
          }
     if(is.infinite(Mo0[[1]])) stop("The posterior is infinite!")
     if(is.nan(Mo0[[1]])) stop("The posterior is not a number!")
     if(is.na(Mo0[[2]])) stop("The deviance is a missing value!")
     if(is.infinite(Mo0[[2]])) stop("The deviance is infinite!")
     if(is.nan(Mo0[[2]])) stop("The deviance is not a number!")
     if(any(is.na(Mo0[[3]])))
          stop("Monitored variable(s) have a missing value!")
     if(any(is.infinite(Mo0[[3]])))
          stop("Monitored variable(s) have an infinite value!")
     if(any(is.nan(Mo0[[3]])))
          stop("Monitored variable(s) include a value that is not a number!")
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
          Fit.LA <- LaplaceApproximation(Model, Initial.Values, Data,
               Method="Rprop", sir=FALSE)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) *
               Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary1[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n")
          cat("for Laplace's Demon, and the posterior modes are now the initial\n")
          cat("values for Laplace's Demon.\n\n")}
     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- matrix(Mo0[[2]], 1, 1)
     Mon <- matrix(Mo0[[3]], 1, length(Mo0[[3]]))
     LIV <- length(Initial.Values)
     post <- matrix(0, Iterations, LIV)
     thinned <- Initial.Values
     post[1,] <- prop <- Initial.Values
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(is.matrix(Covar)) {tuning <- sqrt(diag(Covar)); VarCov <- Covar}
     else {
          tuning <- rep(ScaleF, LIV)
          VarCov <- matrix(0, LIV, LIV)
          diag(VarCov) <- ScaleF / LIV}
     Iden.Mat <- diag(LIV)
     DiagCovar <- matrix(diag(VarCov), 1, LIV)
     ############################  Begin MCMC  ############################
     cat("Algorithm:", Algorithm, "\n")
     cat("\nLaplace's Demon is beginning to update...\n")
     if(Algorithm == "Adaptive Metropolis") {
          mcmc.out <- AM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Adaptive-Mixture Metropolis") {
          mcmc.out <- AMM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, w)}
     else if(Algorithm == "Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- AMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Delayed Rejection Adaptive Metropolis") {
          mcmc.out <- DRAM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Delayed Rejection Metropolis") {
          mcmc.out <- DRM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- MWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Robust Adaptive Metropolis") {
          mcmc.out <- RAM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, alpha.star, Dist, gamma)}
     else if(Algorithm == "Random-Walk Metropolis") {
          mcmc.out <- RWM(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- SAMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, parm.names=Data$parm.names, Dyn)}
     else if(Algorithm == "Sequential Metropolis-within-Gibbs") {
          mcmc.out <- SMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, parm.names=Data$parm.names, Dyn)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- twalk(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, n1=n1, at=at, aw=aw)}
     else if(Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- USAMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, parm.names=Data$parm.names, Dyn, Fit, Begin)}
     else if(Algorithm == "Updating Sequential Metropolis-within-Gibbs") {
          mcmc.out <- USMWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, parm.names=Data$parm.names, Dyn, Fit, Begin)}
     else stop("The algorithm is unrecognized.")
     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     remove(mcmc.out)
     colnames(DiagCovar) <- Data$parm.names
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n")
     ### Real Values
     thinned <- ifelse(!is.finite(thinned), 0, thinned)
     Dev <- ifelse(!is.finite(Dev), 0, Dev)
     Mon <- ifelse(!is.finite(Mon), 0, Mon)
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
          if(class(test) != "try-error") geweke[i,] <- as.vector(test)}
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
          ESS6 <- ESS(Mon[BurnIn:nrow(thinned),])}
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data$parm.names,
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- apply(thinned, 2, sd)
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {Summ1[i,3] <- MCSE(thinned[,i])}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     Deviance[3] <- MCSE(as.vector(Dev))
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          Monitor[3] <- MCSE(Mon[,j])
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.5,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data$mon.names[j]}
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
               Summ2[i,3] <- MCSE(thinned2[,i])}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          Deviance[3] <- MCSE(as.vector(Dev2))
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
               Monitor[3] <- MCSE(as.vector(Mon2[,j]))
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
     if(ncol(Mon) == length(Data$mon.names))
          colnames(Mon) <- Data$mon.names
     if(ncol(thinned) == length(Data$parm.names)) {
          colnames(thinned) <- Data$parm.names}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(({Algorithm == "Delayed Rejection Metropolis"} |
          {Algorithm == "Metropolis-within-Gibbs"} | 
          {Algorithm == "Random-Walk Metropolis"} |
          {Algorithm == "Sequential Metropolis-within-Gibbs"}) &
          {BurnIn < nrow(thinned)}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(theta=thinned2,
               LL=as.vector(Dev2)*(-1/2),
               method="NSIS")}
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,7),
          Adaptive=Adaptive,
          Algorithm=Algorithm,
          Call=LDcall,
          Covar=VarCov,
          CovarDHis=DiagCovar,
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
AM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[[1]]) &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVN.test)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- round(runif(1,0.5,{LIV+0.49}))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
AMM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, w)
     {
     obs.sum <- matrix(0, LIV, 1)
     obs.scatter <- matrix(0, LIV, LIV)
     prop.R <- NULL
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Store Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values from a mixture
          if(is.null(prop.R) || runif(1) < w) {
               prop <- rnorm(LIV, post[iter,], tuning)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Non-Adaptive Component\n")}
          else {
               prop <- t(prop.R) %*% rnorm(LIV) + post[iter,]
               if(iter %% Status == 0) 
                    cat(",   Proposal: Adaptive Component\n")}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + post[iter,]
          obs.scatter <- obs.scatter + tcrossprod(post[iter,])
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               VarCov <- obs.scatter/iter - tcrossprod(obs.sum/iter)
               diag(VarCov) <- diag(VarCov) + 1e-05
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               prop.R <- try(ScaleF * chol(VarCov), silent=TRUE)
               if(!is.matrix(prop.R)) prop.R <- NULL}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
AMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               log.tuning <- ifelse(Acceptance.Rate > 0.44,
                    log.tuning + size, log.tuning - size)
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     VarCov <- cov(post)
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
DRAM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[[1]]) &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVN.test)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- round(runif(1,0.5,{LIV+0.49}))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(is.numeric(MVN.test[[1]]) &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVN.test)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- round(runif(1,0.5,{LIV+0.49}))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo12 <- Model(prop, Data)
               if(!is.finite(Mo12[[1]])) Mo12 <- Mo0
               if(!is.finite(Mo12[[2]])) Mo12 <- Mo0
               if(any(!is.finite(Mo12[[3]]))) Mo12 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[[1]] - Mo12[[1]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo12[[1]] + log.alpha.comp  -
                    {Mo0[[1]] + log(1 - exp(Mo1[[1]] - Mo0[[1]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo12
                    post[iter,] <- Mo12[[5]]
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         Dev[nrow(Dev),] <- Mo1[[2]]
                         Mon[nrow(Mon),] <- Mo1[[3]]}
                    }
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * cov(post[1:iter,])} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
DRM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[[1]]) &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVN.test)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- round(runif(1,0.5,{LIV+0.49}))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- prop
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVN.rand <- rnorm(LIV, 0, 1)
               MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(is.numeric(MVN.test[[1]]) &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVN.test)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- round(runif(1,0.5,{LIV+0.49}))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo12 <- Model(prop, Data)
               if(!is.finite(Mo12[[1]])) Mo12 <- Mo0
               if(!is.finite(Mo12[[2]])) Mo12 <- Mo0
               if(any(!is.finite(Mo12[[3]]))) Mo12 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[[1]] - Mo12[[1]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo12[[1]] + log.alpha.comp  -
                    {Mo0[[1]] + log(1 - exp(Mo1[[1]] - Mo0[[1]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo12
                    post[iter,] <- Mo12[[5]]
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         Dev[nrow(Dev),] <- Mo1[[2]]
                         Mon[nrow(Mon),] <- Mo1[[3]]}
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
MWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Random-Scan Componentwise Estimation
          for (j in sample(LIV)) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
RAM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, alpha.star, Dist, gamma)
     {
     if(!is.symmetric.matrix(VarCov)) {
          cat("\nAsymmetric VarCov, correcting now...\n")
          VarCov <- as.symmetric.matrix(VarCov)}
     if(!is.positive.definite(VarCov)) {
          cat("\nNon-Positive-Definite VarCov, correcting now...\n")
          VarCov <- as.positive.definite(VarCov)}
     S.z <- try(S.z <- t(chol(VarCov)), silent=TRUE)
     if(is.numeric(S.z[[1]])) S <- S.z
     else S <- Iden.Mat
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Store Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose New Parameter Values
          if(Dist == "t") U <- qt(runif(LIV), df = 5, lower.tail = TRUE)
          else U <- rnorm(LIV)
          prop <- post[iter,] + S %*% U
          if(iter %% Status == 0)
               cat(",   Proposal: Multivariate\n")
          ### Log-Posterior
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}}
          ### Adaptation
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               eta <- min(1, LIV*iter^(-gamma))
               VarCov.test <- S %*% (Iden.Mat +
                    eta*(exp(log.alpha) - alpha.star) *
                    U %*% t(U) / sqrt(sum(U^2))) %*% t(S)
               if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
                    !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
               if(!is.symmetric.matrix(VarCov.test))
                    VarCov.test <- as.symmetric.matrix(VarCov.test)
               if(is.positive.definite(VarCov.test)) {
                    S.z <- try(S.z <- t(chol(VarCov)), silent=TRUE)
                    if(is.numeric(S.z[[1]])) {
                         VarCov <- VarCov.test
                         S <- S.z}}
               DiagCovar <- rbind(DiagCovar, diag(VarCov))}
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
RWM <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov)
     {
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter, sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          MVN.rand <- rnorm(LIV, 0, 1)
          MVN.test <- try(MVNz <- matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(is.numeric(MVN.test[[1]]) &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVN.test)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- round(runif(1,0.5,{LIV+0.49}))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- Model(prop, Data)
          if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
          if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
          if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[[1]] - Mo0[[1]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
SAMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, parm.names, Dyn)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          if(length(staticparms) > 1) staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               log.tuning <- ifelse(Acceptance.Rate > 0.44,
                    log.tuning + size, log.tuning - size)
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     VarCov <- cov(post)
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
SMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, parm.names, Dyn)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          if(length(staticparms) > 1) staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
twalk <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, n1, at, aw)
     {
     Obj <- function(x) {return(Model(x, Data)[[1]]*-1)}
     Supp <- function(x) {return(all.equal(x, Model(x, Data)[[5]]))}
     IntProd <- function(x) {return(sum(x*x))}
     DotProd <- function(x, y) {return(sum(x*y))}
     Simh1 <- function(dim, pphi, x, xp, beta)
          {
          phi <- runif(dim) < pphi
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi)))
          }
     Simfbeta <- function(at)
          {
          if(runif(1) < (at-1)/(2*at))
               return(exp(1/(at + 1)*log(runif(1))))
          else
               return(exp(1/(1 - at)*log(runif(1))))
          }
     Simh2 <- function(dim, pphi, aw, x, xp)
          {
          u <- runif(dim)
          phi <- runif(dim) < pphi
          z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
          z <- z*phi
          return(list(rt=x + (x - xp)*z, nphi=sum(phi)))
          }
     Simh3 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))
          x + sigma*rnorm(dim)*phi
          return(list(rt=x + sigma*rnorm(dim)*phi, nphi=sum(phi),
               sigma=sigma))
          }
     G3U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd(h - xp)/(sigma^2))
          else
               return(0) 
          }
     Simh4 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))/3
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + sigma*rnorm(1))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi), sigma=sigma))
          }
     G4U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd((h - x))/(sigma^2))
          else
               return(0)
          }
     OneMove <- function(dim, Obj, Supp, x, U, xp, Up, at=at, aw=aw,
          pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918)
          {
          ker <- runif(1) ### Choose a kernel
          if(ker < F1) {
               ### Kernel h1: traverse	
               dir <- runif(1) 
               funh <- 1
               if((0 <= dir) && (dir < 0.5)) {	
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, xp, x, beta)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    if(Supp(yp)) {
                         propUp <- Obj(yp)
                         ### The proposal is symmetric
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               if((0.5 <= dir) && (dir < 1.0)) {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, x, xp, beta)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    if(Supp(y)) {
                         propU <- Obj(y)
                         ### The proposal is symmetric
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if((F1 <= ker) && (ker < F2)) {
               ### Kernel h2: walk
               dir <- runif(1)
               funh <- 2
               if((0 <= dir) && (dir < 0.5)) {
                    ### x as pivot
                    tmp <- Simh2(dim, pphi, aw, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    if((Supp(yp)) && (all(abs(yp - y) > 0))) {
                         propUp <- Obj(yp)
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               if((0.5 <= dir) && (dir < 1.0)) {
                    ### xp as pivot
                    tmp <- Simh2(dim, pphi, aw, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    if((Supp(y)) && (all(abs(yp - y) > 0))) {
                         propU <- Obj(y)
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if((F2 <= ker) && (ker < F3)) {
               ### Kernel h3: blow
               dir <- runif(1)
               funh <- 3
               if((0 <= dir) && (dir < 0.5)) {
                    ### x as pivot
                    tmp <- Simh3(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    if((Supp(yp)) && all(yp != x)) {
                         propUp <- Obj(yp)
                         W1 <- G3U(nphi, sigma,  yp, xp,  x)
                         W2 <- G3U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               if((0.5 <= dir) && (dir < 1.0)) {
                    ### xp as pivot
                    tmp <- Simh3(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    if((Supp(y)) && all(y != xp)) {
                         propU <- Obj(y)
                         W1 <- G3U(nphi, sigma,   y,  x, xp)
                         W2 <- G3U(nphi, sigma,   x,  y, xp)
                         A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(F3 <= ker) {
               ## Kernel h4: hop
               dir <- runif(1)
               funh <- 4
               if((0 <= dir) && (dir < 0.5)) {
                    ### x as pivot
                    tmp <- Simh4(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    if((Supp(yp)) && all(yp != x)) {
                         propUp <- Obj(yp)
                         W1 <- G4U(nphi, sigma,  yp, xp,  x)
                         W2 <- G4U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               if((0.5 <= dir) && (dir < 1.0)) {
                    ### xp as pivot
                    tmp <- Simh4(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    if((Supp(y)) && all(y != xp)) {
                         propU <- Obj(y)
                         W1 <- G4U(nphi, sigma,   y,  x, xp)
                         W2 <- G4U(nphi, sigma,   x,  y, xp)
                         A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(is.nan(A)) {
               #### Debugging
               cat("twalk: Error in evaluating the objective:")
               cat( funh, "DU=", (U - propU), "DUp=", (Up - propUp),
                    "DW=", (W1 - W2), "A=", A, "\n", y, "\n", yp, "\n")
               cat("U=", U, "propU=", propU, "Up=", Up, "propUp", propUp)}
          return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A,
               funh=funh, nphi=nphi))
          }
     Runtwalk <- function(Iterations, dim, Obj, Supp, x0, xp0, 
          pphi, at, aw, F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model,
          Data, Status, Thinning, Acceptance, Dev, Mon, Mo0, post, thinned)
          {
          x <- x0
          xp <- xp0
          if(Supp(x) && Supp(xp)) {
               U <- Obj(x)
               Up <- Obj(xp)
               }
          else {
               cat(paste("Initial values out of support,\n  x=", x,"\n xp=", xp))
               Iterations <- 0}
          if(any(abs(x0 - xp0) <= 0)) {
               cat("\nBoth sets of initial values are not unique.")
               Iterations <- 0}
          Acceptance <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) {cat("Iteration: ", iter,
                    ",   Proposal: Subset Multivariate\n", sep="")}
               ### Store Current Posterior
               if(iter > 1) post[iter,] <- post[iter-1,]
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    thinned <- rbind(thinned, post[iter,])
                    Dev <- rbind(Dev, Mo0[[2]])
                    Mon <- rbind(Mon, Mo0[[3]])}
               ### Propose New Parameter Values
               move <- OneMove(dim=dim, Obj=Obj, Supp=Supp, x, U, xp, Up,
                    at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3)
               ### Log-Posterior of the proposed state
               Mo1 <- Model(x, Data)
               ### Accept/Reject
               if(runif(1) < move$A) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo1[[5]]
                    Acceptance <- Acceptance + 1 #move$nphi/dim
                    if(iter %% Thinning == 0) {
                         Dev[nrow(Dev),] <- Mo1[[2]]
                         Mon[nrow(Mon),] <- Mo1[[3]]}
                    x <- move$y
                    U <- move$propU
                    xp <- move$yp
                    Up <- move$propUp
                    }
               }
          out <- list(Acceptance=Acceptance,
               Dev=Dev,
               DiagCovar=DiagCovar,
               Mon=Mon,
               thinned=thinned,
               VarCov=VarCov)
          return(out)
          }
     out <- Runtwalk(Iterations=Iterations, dim=LIV, Obj=Obj, Supp=Supp,
          x0=post[1,], xp0=GIV(Model, Data), pphi=n1/LIV, at=6, aw=1.5,
          Model=Model, Data=Data, Status=Status, Thinning=Thinning,
          Acceptance=Acceptance, Dev=Dev, Mon=Mon, Mo0=Mo0, post=post,
          thinned=thinned)
     ### Output
     return(out)
     }
USAMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, parm.names, Dyn, Fit,
     Begin)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.")
     ivs <- post[1,]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               log.tuning <- ifelse(Acceptance.Rate > 0.44,
                    log.tuning + size, log.tuning - size)
               tuning <- exp(log.tuning)
               DiagCovar <- rbind(DiagCovar, tuning)}
          }
     VarCov <- cov(post)
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
USMWG <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, parm.names, Dyn, Fit,
     Begin)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.")
     ivs <- post[1,]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) {cat("Iteration: ", iter,
               ",   Proposal: Componentwise\n", sep="")}
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          if(ncol(Dyn) > 1) dynsample <- as.vector(apply(Dyn,1,sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- Model(prop, Data)
               if(!is.finite(Mo1[[1]])) Mo1 <- Mo0
               if(!is.finite(Mo1[[2]])) Mo1 <- Mo0
               if(any(!is.finite(Mo1[[3]]))) Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[[1]] - Mo0[[1]])
               post[iter,j] <- Mo1[[5]][j]*(u == 1) + post[iter,j]*(u == 0)
               Mo0[[1]] <- Mo1[[1]]*(u == 1) + Mo0[[1]]*(u == 0)
               Mo0[[2]] <- Mo1[[2]]*(u == 1) + Mo0[[2]]*(u == 0)
               Mo0[[3]] <- Mo1[[3]]*(u == 1) + Mo0[[3]]*(u == 0)
               Acceptance[j] <- Acceptance[j] + u}
          if(iter %% Thinning == 0) {
               Dev[nrow(Dev),] <- Mo1[[2]]
               Mon[nrow(Mon),] <- Mo1[[3]]}
          }
     VarCov <- cov(post)
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
          Dev=Dev,
          DiagCovar=DiagCovar,
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }

#End
