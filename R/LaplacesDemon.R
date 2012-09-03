###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=100000, Status=1000, Thinning=100, Algorithm="RWM",
     Specs=NULL, ...)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="")
     time1 <- proc.time()
     LDcall <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n")
     if(missing(Model)) stop("A function must be entered for Model.")
     if(!is.function(Model)) stop("Model must be a function.")
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
          cat("WARNING: The length of Initial Values differed from",
               "Data$parm.names.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n")
          Initial.Values <- rep(0, length(Data$parm.names))}
     Iterations <- round(abs(Iterations))
     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="")}
     Status <- round(abs(Status))
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="")}
     Thinning <- round(abs(Thinning))
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="")}
     if(Algorithm %in% c("AHMC","AM","AMM","AMWG","DRAM","DRM",
          "Experimental","HMC","HMCDA","MWG","NUTS","RAM","RWM","SAMWG",
          "SMWG","THMC","twalk","USAMWG","USMWG")) {
          if(Algorithm == "AHMC") {
               Algorithm <- "Adaptive Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- as.vector(abs(Specs[[1]]))
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 1) {
                    cat("\nL has been increased to its minimum: 1.\n")
                    L <- 1}
               Adaptive <- Iterations + 1
               DR <- 1
               if(Specs[["Periodicity"]] == Specs[[3]]) Periodicity <- Specs[[3]]}
          if(Algorithm == "AM") {
               Algorithm <- "Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}
               if(Specs[["w"]] == Specs[[3]]) {
                    w <- Specs[[3]]
                    if(w <= 0 || w >= 1) w <- 0.05}
               else {
                    w <- 0.05
                    cat("\nw was misspecified and changed to 0.05.\n")}}
          if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 1) stop("The Specs argument is incorrect.")
               Adaptive <- 2
               DR <- 0
               if(Specs[["Periodicity"]] == Specs[[1]]) Periodicity <- Specs[[1]]
               else {
                    Periodicity <- 50
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}}
          if(Algorithm == "DRAM") {
               Algorithm <- "Delayed Rejection Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(Specs[["Adaptive"]] == Specs[[1]]) Adaptive <- Specs[[1]]
               else {
                    Adaptive <- 20
                    cat("\nAdaptive was misspecified and changed to 20.\n")}
               DR <- 1
               if(Specs[["Periodicity"]] == Specs[[2]]) Periodicity <- Specs[[2]]
               else {
                    Periodicity <- 100
                    cat("\nPeriodicity was misspecified and changed to 100.\n")}}
          if(Algorithm == "DRM") {
               Algorithm <- "Delayed Rejection Metropolis"
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "Experimental") {
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               Adaptive <- Iterations + 1
               DR <- 0
               Periodicity <- Iterations + 1}
          if(Algorithm == "HMC") {
               Algorithm <- "Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 2) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- abs(Specs[[1]])
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 1) {
                    cat("\nL has been increased to its minimum: 1.\n")
                    L <- 1}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "HMCDA") {
               Algorithm <- "Hamiltonian Monte Carlo with Dual-Averaging"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 5) stop("The Specs argument is incorrect.")
               if(Specs[["A"]] == Specs[[1]])
                    A <- min(round(abs(Specs[[1]])), Iterations)
               if(Specs[["delta"]] == Specs[[2]])
                    delta <- max(min(abs(Specs[[2]]), 1), 1/Iterations)
               if(is.null(Specs[[3]]) |
                    {all(Specs[["epsilon"]] == Specs[[3]])}) {
                    if(is.null(Specs[[3]])) epsilon <- NULL
                    else epsilon <- abs(Specs[[3]][1])}
               if(Specs[["Lmax"]] == Specs[[4]]) 
                    Lmax <- abs(Specs[[4]])
               if(Specs[["lambda"]] == Specs[[5]]) {
                    lambda <- abs(Specs[[5]])
                    if(!is.null(epsilon))
                         if(lambda < epsilon) lambda <- epsilon}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "NUTS") {
               Algorithm <- "No-U-Turn Sampler"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(Specs[["A"]] == Specs[[1]])
                    A <- max(min(round(abs(Specs[[1]])), Iterations),1)
               if(Specs[["delta"]] == Specs[[2]])
                    delta <- max(min(abs(Specs[[2]]), 1), 1/Iterations)
               if(is.null(Specs[[3]]) |
                    {all(Specs[["epsilon"]] == Specs[[3]])}) {
                    if(is.null(Specs[[3]])) epsilon <- NULL
                    else epsilon <- abs(Specs[[3]][1])}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- Iterations + 1}
          if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               if(Specs[["alpha.star"]] == Specs[[1]]) {
                    alpha.star <- Specs[[1]]
                    if(alpha.star <= 0 || alpha.star >= 1) {
                         cat("\nalpha.star not in (0,1). Changed to 0.234.\n")
                         alpha.star <- 0.234}
                    if({length(alpha.star) != 1} &
                         {length(alpha.star) != length(Initial.Values)}) {
                         cat("\nLength of alpha.star is wrong. Changed to 1.\n")
                         alpha.star <- alpha.star[1]}}
               else {
                    alpha.star <- 0.234
                    cat("\nalpha.star was misspecified and changed to 0.234.\n")}
               Adaptive <- 2
               DR <- 0
               if(Specs[["Dist"]] == Specs[[2]]) {
                    Dist <- Specs[[2]]
                    if(Dist != "t" & Dist != "N") {
                         cat("\nDist was not t or N, and changed to N.\n")
                         Dist <- "N"}}
               else {
                    Dist <- "N"
                    cat("\nDist was not t or N, and changed to N.\n")
                    }
               if(Specs[["gamma"]] == Specs[[3]]) {
                    gamma <- Specs[[3]]
                    if(gamma <= 0.5 || gamma > 1) {
                         cat("\ngamma not in (0.5,1]. Changed to 0.66.\n")
                         gamma <- 0.66}}
               else {
                    gamma <- 0.66
                    cat("\ngamma was misspecified and changed to 0.66.\n")}
               if(Specs[["Periodicity"]] == Specs[[4]]) Periodicity <- Specs[[4]]
               else {
                    Periodicity <- 10
                    cat("\nPeriodicity was misspecified and changed to 10.\n")}}
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
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}}
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
          if(Algorithm == "THMC") {
               Algorithm <- "Tempered Hamiltonian Monte Carlo"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 3) stop("The Specs argument is incorrect.")
               if(all(Specs[["epsilon"]] == Specs[[1]])) {
                    epsilon <- as.vector(abs(Specs[[1]]))
                    if(length(epsilon) != length(Initial.Values))
                         cat("\nLength of epsilon is incorrect.\n")
                         epsilon <- rep(epsilon[1], length(Initial.Values))}
               if(Specs[["L"]] == Specs[[2]]) L <- abs(round(Specs[[2]]))
               if(L < 2) {
                    cat("\nL has been increased to its minimum: 2.\n")
                    L <- 2}
               Adaptive <- Iterations + 1
               DR <- 1
               Periodicity <- 1
               if(Specs[["Temperature"]] == Specs[[3]])
                    Temperature <- Specs[[3]]
                    if(Temperature <= 0) {
                         cat("\nTemperature is incorrect, changed to 1.\n")
                         Temperature <- 1}}
          if(Algorithm == "twalk") {
               Algorithm <- "t-walk"
               if(missing(Specs)) stop("The Specs argument is required.")
               if(length(Specs) != 4) stop("The Specs argument is incorrect.")
               Adaptive <- Iterations + 1
               DR <- 0
               if(is.null(Specs[[1]]) | all(Specs[["SIV"]] == Specs[[1]])) {
                    if(is.null(Specs[[1]])) {
                         cat("\nGenerating SIV...\n")
                         if(!is.null(Data$PGF))
                              SIV <- GIV(Model, Data, PGF=TRUE)
                         else SIV <- GIV(Model, Data)}
                    else SIV <- Specs[[1]]
                    if(!identical(length(SIV), length(Initial.Values))) {
                         cat("\nGenerating SIV due to length mismatch.\n")
                         if(!is.null(Data$PGF))
                              SIV <- GIV(Model, Data, PGF=TRUE)
                         else SIV <- GIV(Model, Data)}}
               else {
                    cat("\nSIV was misspecified. Generating values now.\n")
                    if(!is.null(Data$PGF))
                         SIV <- GIV(Model, Data, PGF=TRUE)
                    else SIV <- GIV(Model, Data)}
               Mo2 <- Model(SIV, Data)
               if(!is.finite(Mo2[[1]]))
                    stop("SIV results in a non-finite posterior.")
               if(!is.finite(Mo2[[2]]))
                    stop("SIV results in a non-finite deviance.")
               SIV <- Mo2[[5]]
               rm(Mo2)
               if(Specs[["n1"]] == Specs[[2]]) {
                    n1 <- Specs[[2]]
                    if(n1 < 1) {
                         cat("\nn1 must be at least 1. Changed to 4.\n")
                         n1 <- 4}}
               else {
                    n1 <- 4
                    cat("\nn1 was misspecified and changed to 4.\n")}
               if(Specs[["at"]] == Specs[[3]]) {
                    at <- Specs[[3]]
                    if(at <= 0) {
                         cat("\nat must be positive. Changed to 6.\n")
                         at <- 6}}
               else {
                    at <- 6
                    cat("\nat was misspecified and changed to 6.\n")}
               if(Specs[["aw"]] == Specs[[4]]) {
                    aw <- Specs[[4]]
                    if(aw <= 0) {
                         cat("\naw must be positive. Changed to 1.5.\n")
                         at <- 1.5}}
               else {
                    aw <- 1.5
                    cat("\naw was misspecified and changed to 1.5.\n")}
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
                    cat("\nPeriodicity was misspecified and changed to 50.\n")}
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
     if(!is.list(Mo0)) stop("Model must return a list.")
     if(length(Mo0) != 5) stop("Model must return five components.")
     if(length(Mo0[[1]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(Mo0[[3]]), length(Data$mon.names)))
          stop("Length of mon.names differs from length of monitors.")
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
          if(!is.null(Data$PGF))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
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
     if(Algorithm == "t-walk") {
          Mo0 <- Model(Initial.Values, Data)
          if(any(Mo0[[5]] == SIV))
              stop("Initial.Values and SIV not unique after model update.")}
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
     if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- AHMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, epsilon, L)}
     else if(Algorithm == "Adaptive Metropolis") {
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
     else if(Algorithm == "Experimental") {
          stop("Experimental function not found.")}
     else if(Algorithm == "Hamiltonian Monte Carlo") {
          mcmc.out <- HMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, epsilon, L)}
     else if(Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          mcmc.out <- HMCDA(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, A, delta, epsilon, Lmax, lambda)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- MWG(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov)}
     else if(Algorithm == "No-U-Turn Sampler") {
          mcmc.out <- NUTS(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, A, delta, epsilon)}
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
     else if(Algorithm == "Tempered Hamiltonian Monte Carlo") {
          mcmc.out <- THMC(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, epsilon, L, Temperature)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- twalk(Model, Data, Adaptive, DR, Iterations,
               Periodicity, Status, Thinning, Acceptance, Dev, DiagCovar,
               Iden.Mat, LIV, Mon, Mo0, post, ScaleF, thinned, tuning,
               VarCov, SIV=SIV, n1=n1, at=at, aw=aw)}
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
     colnames(VarCov) <- rownames(VarCov) <- Data$parm.names
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
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)))
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
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(Mon[,j]), silent=TRUE)
          if(inherits(temp, "try-error")) 
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
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
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i],
                    method="sample.variance")}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- MCSE(as.vector(Dev2))
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
          {Algorithm == "Hamiltonian Monte Carlo"} | 
          {Algorithm == "Metropolis-within-Gibbs"} |
          {Algorithm == "No-U-Turn Sampler"} | 
          {Algorithm == "Random-Walk Metropolis"} |
          {Algorithm == "Sequential Metropolis-within-Gibbs"} |
          {Algorithm == "Tempered Hamiltonian Monte Carlo"} | 
          {Algorithm == "t-walk"}) &
          {BurnIn < nrow(thinned)}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(theta=thinned2,
               LL=as.vector(Dev2)*(-1/2), method="NSIS")}
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
AHMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, epsilon, L)
     {
     DiagCovar[1,] <- epsilon
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[[5]] == Mo1[[5]])) {
                    nomove <- which(Mo0.1[[5]] == Mo1[[5]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[[5]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[[1]] + kinetic0
          H1 <- -Mo1[[1]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Adaptation
          if({iter > 10} & {iter %% Periodicity == 0}) {
               acceptances <- apply(post[(iter-9):iter,], 2, function(x)
                    {length(unique(x))})
               eps.num <- which(acceptances <= 1)
               epsilon[eps.num] <- epsilon[eps.num] * 0.8
               eps.num <- which(acceptances > 7)
               epsilon[eps.num] <- epsilon[eps.num] * 1.2
               DiagCovar <- rbind(DiagCovar, epsilon)}
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
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
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
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
               MVNz <- try(matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
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
               MVNz <- try(matrix(MVN.rand,1,LIV) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
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
HMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, epsilon, L)
     {
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[[5]] == Mo1[[5]])) {
                    nomove <- which(Mo0.1[[5]] == Mo1[[5]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[[5]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[[1]] + kinetic0
          H1 <- -Mo1[[1]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
HMCDA <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, A, delta, epsilon, Lmax,
     lambda)
     {
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
          thetaprime <- Mo1[[5]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data)
          {
          cat("\nFinding a reasonable initial value for epsilon...")
          epsilon <- 1
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
          acceptprob <- exp(leap$Mo1[[1]] - Mo0[[1]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
               acceptprob <- exp(leap$Mo1[[1]] - Mo0[[1]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="")
          return(epsilon)
          }
     gr0 <- partial(Model, post[1,], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(post[1,], gr0, Mo0, Model,
               Data)
     DiagCovar[1,] <- epsilon
     L <- max(1, round(lambda / epsilon))
     L <- min(L, Lmax)
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Begin HMCDA
     for (iter in 1:Iterations) {
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum1 <- momentum0 <- runif(LIV)
          joint <- Mo0[[1]] - 0.5 * as.vector(momentum0 %*% momentum0)
          L <- max(1, round(lambda / epsilon))
          L <- min(L, Lmax)
          gr1 <- gr0
          Mo0.1 <- Mo0
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate,   L: ", L, "\n", sep="")
          ### Leapfrog Function
          for (l in 1:L) {
               momentum1 <- momentum1 + 0.5 * epsilon * gr1
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[[5]] == Mo1[[5]])) {
                    nomove <- which(Mo0.1[[5]] == Mo1[[5]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[[5]]
               gr1 <- partial(Model, prop, Data)
               momentum1 <- momentum1 + epsilon * gr1}
          ### Accept/Reject
          alpha <- min(1,
               exp(prop - 0.5 * as.vector(momentum1 %*% momentum1) - joint))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          ### Adaptation
          if(iter > 1) {
               eta <- 1 / (iter - 1 + t0)
               Hbar <- (1 - eta) * Hbar + eta * (delta - alpha)
               if(iter <= A) {
                    epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
                    eta <- (iter-1)^-kappa
                    epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                         eta * log(epsilon))
                    DiagCovar <- rbind(DiagCovar, epsilon)}
               else epsilon <- epsilonbar}
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
NUTS <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, A, delta, epsilon)
     {
     leapfrog <- function(theta, r, grad, epsilon, Model, Data)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- Model(thetaprime, Data)
          thetaprime <- Mo1[[5]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime,
               rprime=rprime,
               gradprime=gradprime,
               Mo1=Mo1)
          return(out)
          }
     stop.criterion <- function(thetaminus, thetaplus, rminus, rplus)
          {
          thetavec <- thetaplus - thetaminus
          criterion <- (thetavec %*% rminus >= 0) &&
               (thetavec %*% rplus >= 0)
          return(criterion)
          }
     build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0)
          {
          if(j == 0) {
               ### Base case: Take a single leapfrog step in direction v
               leap <- leapfrog(theta, r, grad, v*epsilon, Model, Data)
               rprime <- leap$rprime
               thetaprime <- leap$thetaprime
               Mo1 <- leap$Mo1
               gradprime <- leap$gradprime
               joint <- Mo1[[1]] - 0.5 * as.vector(rprime %*% rprime)
               ### Is the new point in the slice?
               nprime <- logu < joint
               ### Is the simulation wildly inaccurate?
               sprime <- logu - 1000 < joint
               # Set the return values---minus=plus for all things here,
               # since the "tree" is of depth 0.
               thetaminus <- thetaprime
               thetaplus <- thetaprime
               rminus <- rprime
               rplus <- rprime
               gradminus <- gradprime
               gradplus <- gradprime
               ### Compute the acceptance probability
               alphaprime <- min(1, exp(Mo1[[1]] - 0.5 *
                    as.vector(rprime %*% rprime) - joint0))
               nalphaprime <- 1}
          else {
               # Recursion: Implicitly build the height j-1 left and
               # right subtrees
               tree <- build.tree(theta, r, grad, logu, v, j-1, epsilon,
                    joint0)
               thetaminus <- tree$thetaminus
               rminus <- tree$rminus
               gradminus <- tree$gradminus
               thetaplus <- tree$thetaplus
               rplus <- tree$rplus
               gradplus <- tree$gradplus
               thetaprime <- tree$thetaprime
               gradprime <- tree$gradprime
               Mo1 <- tree$Mo1
               nprime <- tree$nprime
               sprime <- tree$sprime
               alphaprime <- tree$alphaprime
               nalphaprime <- tree$nalphaprime
               ### If the first subtree stopping criterion is met, then stop
               if(sprime == 1) {
                    if(v == -1) {
                         tree <- build.tree(thetaminus, rminus, gradminus,
                              logu, v, j-1, epsilon, joint0)
                         thetaminus <- tree$thetaminus
                         rminus <- tree$rminus
                         gradminus <- tree$gradminus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    else {
                         tree <- build.tree(thetaplus, rplus, gradplus,
                              logu, v, j-1, epsilon, joint0)
                         thetaplus <- tree$thetaplus
                         rplus <- tree$rplus
                         gradplus <- tree$gradplus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    ### Choose a subtree to propagate a sample up from
                    temp <- nprime2 / (nprime + nprime2)
                    if(!is.finite(temp)) temp <- 0
                    if(runif(1) < temp) {
                         thetaprime <- thetaprime2
                         gradprime <- gradprime2
                         Mo1 <- Mo12}
                    ### Update the number of valid points
                    nprime <- nprime + nprime2
                    ### Update the stopping criterion
                    sprime <- sprime && sprime2 &&
                         stop.criterion(thetaminus, thetaplus, rminus,
                              rplus)
                    ### Update acceptance probability statistics
                    alphaprime <- alphaprime + alphaprime2
                    nalphaprime <- nalphaprime + nalphaprime2}}
          out <- list(thetaminus=thetaminus,
               rminus=rminus,
               gradminus=gradminus,
               thetaplus=thetaplus,
               rplus=rplus,
               gradplus=gradplus,
               thetaprime=thetaprime,
               gradprime=gradprime,
               Mo1=Mo1,
               nprime=nprime,
               sprime=sprime,
               alphaprime=alphaprime,
               nalphaprime=nalphaprime)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data)
          {
          cat("\nFinding a reasonable initial value for epsilon...")
          epsilon <- 1
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
          acceptprob <- exp(leap$Mo1[[1]] - Mo0[[1]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data)
               acceptprob <- exp(leap$Mo1[[1]] - Mo0[[1]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="")
          return(epsilon)
          }
     Count <- 0
     evals <- 0
     grad <- partial(Model, post[1,], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(post[1,], grad, Mo0, Model,
               Data)
     DiagCovar[1,] <- epsilon
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Begin NUTS
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0 & {iter > A | A >= Iterations}) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          prop <- post[iter,]
          r0 <- runif(LIV) ### r0 is momenta
          ### Joint log-probability of theta and momenta r
          joint <- Mo0[[1]] - 0.5 * as.vector(r0 %*% r0)
          ### Resample u ~ U([0, exp(joint)])
          logu <- joint - rexp(1)
          ### Initialize Tree
          thetaminus <- prop
          thetaplus <- prop
          rminus <- r0
          rplus <- r0
          gradminus <- grad
          gradplus <- grad
          j <- 0 ### Initial height j=0
          n <- 1 ### Initially, the only valid point is the initial point
          s <- 1 ### Loop until s == 0
          while (s == 1) {
               ### Choose a direction: -1=backwards, 1=forwards.
               v <- 2*(runif(1) < 0.5) - 1
               ### Double the size of the tree.
               if(v == -1) {
                    tree <- build.tree(thetaminus, rminus, gradminus,
                         logu, v, j, epsilon, joint)
                    thetaminus <- tree$thetaminus
                    rminus <- tree$rminus
                    gradminus <- tree$gradminus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               else {
                    tree <- build.tree(thetaplus, rplus, gradplus, logu,
                         v, j, epsilon, joint)
                    thetaplus <- tree$thetaplus
                    rplus <- tree$rplus
                    gradplus <- tree$gradplus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               ### Accept/Reject
               Count <- Count + 1
               if((sprime == 1) && (runif(1) < nprime/n)) {
                      post[iter,] <- thetaprime
                      Mo0 <- Mo1
                      grad <- gradprime
                      Acceptance <- Acceptance + 1
                      if(iter %% Thinning == 0 &
                           {iter > A | A >= Iterations}) {
                           #Acceptance <- Acceptance + 1
                           Dev[nrow(Dev),] <- Mo1[[2]]
                           Mon[nrow(Mon),] <- Mo1[[3]]}}
               ### Update number of observed valid points
               n <- n + nprime
               ### Decide if it is time to stop
               s <- sprime &&
                    stop.criterion(thetaminus, thetaplus, rminus, rplus)
               ### Increment depth
               j <- j + 1}
          ### Adaptation of epsilon
          eta <- 1 / (iter - 1 + t0)
          Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
          if(iter <= A) {
               epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
               eta <- (iter-1)^-kappa
               epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                    eta * log(epsilon))
               DiagCovar <- rbind(DiagCovar, epsilon)}
          else epsilon <- epsilonbar
          }
     Acceptance <- round(Acceptance / Count * Iterations)
     ### Output
     out <- list(Acceptance=Acceptance,
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
     S.z <- try(t(chol(VarCov)), silent=TRUE)
     if(!inherits(S.z, "try-error")) S <- S.z
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
          if(Dist == "t") U <- qt(runif(LIV), df=5, lower.tail=TRUE)
          else U <- rnorm(LIV)
          prop <- post[iter,] + S %*% U
          if(iter %% Status == 0)
               cat(",   Proposal: Multivariate\n")
          ### Log-Posterior
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
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
                    S.z <- try(t(chol(VarCov)), silent=TRUE)
                    if(!inherits(S.z, "try-error")) {
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
          MVNz <- try(matrix(MVN.rand,1,LIV) %*% chol(VarCov),
               silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate\n")
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component\n")
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
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
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
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
THMC <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, epsilon, L, Temperature)
     {
     gr <- partial(Model, post[1,], Data)
     sqrt.Temp <- sqrt(Temperature)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0) cat("Iteration: ", iter,
               ",   Proposal: Multivariate\n", sep="")
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               thinned <- rbind(thinned, post[iter,])
               Dev <- rbind(Dev, Mo0[[2]])
               Mon <- rbind(Mon, Mo0[[3]])}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum1 <- momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          Mo0.1 <- Mo0
          for (l in 1:L) {
               if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
               else momentum1 <- momentum1 / sqrt.Temp
               momentum1 <- momentum1 + (epsilon/2) * gr
               prop <- prop + epsilon * momentum1
               Mo1 <- Model(prop, Data)
               if(any(Mo0.1[[5]] == Mo1[[5]])) {
                    nomove <- which(Mo0.1[[5]] == Mo1[[5]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- Model(prop, Data)}
               Mo0.1 <- Mo1
               prop <- Mo1[[5]]
               gr <- partial(Model, prop, Data)
               momentum1 <- momentum1 + (epsilon/2) * gr
               if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
               else momentum1 <- momentum1 * sqrt.Temp}
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[[1]] + kinetic0
          H1 <- -Mo1[[1]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[[5]]
               kinetic0 <- kinetic1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    Dev[nrow(Dev),] <- Mo1[[2]]
                    Mon[nrow(Mon),] <- Mo1[[3]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance,
          Dev=Dev,
          DiagCovar=matrix(epsilon, 1, LIV),
          Mon=Mon,
          thinned=thinned,
          VarCov=VarCov)
     return(out)
     }
twalk <- function(Model, Data, Adaptive, DR, Iterations, Periodicity,
     Status, Thinning, Acceptance, Dev, DiagCovar, Iden.Mat, LIV, Mon,
     Mo0, post, ScaleF, thinned, tuning, VarCov, SIV, n1, at, aw)
     {
     xp0 <- SIV
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
     OneMove <- function(dim, Model, Data, x, U, xp, Up, at=at, aw=aw,
          pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918, Mo0.1, Mo0.2)
          {
          dir <- runif(1) ### Determine which set of points
          ker <- runif(1) ### Choose a kernel
          if(ker < F1) {
               ### Kernel h1: Traverse
               funh <- 1
               if(dir < 0.5) {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, xp, x, beta)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[[1]]) &
                              identical(yp, as.vector(Mo1.2[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propUp <- Mo1.2[[1]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, x, xp, beta)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[[1]]) &
                              identical(y, as.vector(Mo1.1[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propU <- Mo1.1[[1]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F2) {
               ### Kernel h2: Walk
               funh <- 2
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh2(dim, pphi, aw, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[[1]]) &
                              identical(yp, as.vector(Mo1.2[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propUp <- Mo1.2[[1]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh2(dim, pphi, aw, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[[1]]) &
                              identical(y, as.vector(Mo1.1[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propU <- Mo1.1[[1]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F3) {
               ### Kernel h3: Blow
               funh <- 3
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh3(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[[1]]) &
                              identical(yp, as.vector(Mo1.2[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[[1]] * -1
                         W1 <- G3U(nphi, sigma,  yp, xp,  x)
                         W2 <- G3U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh3(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[[1]]) &
                              identical(y, as.vector(Mo1.1[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[[1]] * -1
                         W1 <- G3U(nphi, sigma, y, x, xp)
                         W2 <- G3U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else {
               ## Kernel h4: Hop
               funh <- 4
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh4(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[[1]]) &
                              identical(yp, as.vector(Mo1.2[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[[1]] * -1
                         W1 <- G4U(nphi, sigma, yp, xp, x)
                         W2 <- G4U(nphi, sigma, xp, yp, x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh4(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[[1]]) &
                              identical(y, as.vector(Mo1.1[[5]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[[1]] * -1
                         W1 <- G4U(nphi, sigma, y, x, xp)
                         W2 <- G4U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(check1 & check2 & is.finite(A) & (dir < 0.5))
               Mo0.2 <- Mo1.2
          else if(check1 & check2 & is.finite(A) & (dir >= 0.5))
               Mo0.1 <- Mo1.1
          else if(!is.finite(A)) A <- 0
          #else if(!is.finite(A)) {
          #     #### Debugging
          #     cat("\ntwalk: Error in evaluating the objective:\n")
          #     cat( funh, "DU=", (U - propU), "DUp=", (Up - propUp),
          #          "DW=", (W1 - W2), "A=", A, "\n", y, "\n", yp, "\n")
          #     cat("U=", U, "propU=", propU, "Up=", Up, "propUp", propUp)}
          return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A,
               funh=funh, nphi=nphi, Mo0.1=Mo0.1, Mo0.2=Mo0.2))
          }
     Runtwalk <- function(Iterations, dim, x0, xp0, pphi, at, aw,
          F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model, Data, Status,
          Thinning, Acceptance, Dev, Mon, Mo0, post, thinned, VarCov)
          {
          x <- x0 ### Primary vector of initial values
          xp <- xp0 ### Secondary vector of initial values
          Mo0.1 <- try(Model(x, Data), silent=TRUE)
          Mo0.2 <- try(Model(xp, Data), silent=TRUE)
          if(inherits(Mo0.1, "try-error") | inherits(Mo0.2, "try-error"))
               stop("Error in estimating the log-posterior.")
          if(!is.finite(Mo0.1[[1]]) | !is.finite(Mo0.2[[1]]))
               stop("The log-posterior is non-finite.")
          if(identical(x, as.vector(Mo0.1[[5]])) &
               identical(xp, as.vector(Mo0.2[[5]]))) {
               U <- Mo0.1[[1]] * -1
               Up <- Mo0.2[[1]] * -1}
          else {
               cat("\nInitial values are out of support.")
               cat("\n  Initial.Values=", x)
               cat("\n SIV=", xp)
               stop("Try re-specifying initial values.")}
          if(any(abs(x - xp) <= 0))
               stop("\nBoth vectors of initial values are not unique.")
          Acceptance <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0) {cat("Iteration: ", iter,
                    ",   Proposal: Multivariate Subset\n", sep="")}
               ### Store Current Posterior
               if(iter > 1) post[iter,] <- post[iter-1,]
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    thinned <- rbind(thinned, post[iter,])
                    Dev <- rbind(Dev, Mo0.1[[2]])
                    Mon <- rbind(Mon, Mo0.1[[3]])}
               ### Assign x and xp
               x <- as.vector(Mo0.1[[5]])
               xp <- as.vector(Mo0.2[[5]])
               ### Propose New Parameter Values
               move <- OneMove(dim=dim, Model, Data, x, U, xp, Up,
                    at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3,
                    Mo0.1=Mo0.1, Mo0.2=Mo0.2)
               ### Accept/Reject
               if(runif(1) < move$A) {
                    Mo0.1 <- move$Mo0.1
                    Mo0.2 <- move$Mo0.2
                    post[iter,] <- move$Mo0.1[[5]]
                    Acceptance <- Acceptance + 1 #move$nphi/dim
                    if(iter %% Thinning == 0) {
                         Dev[nrow(Dev),] <- move$Mo0.1[[2]]
                         Mon[nrow(Mon),] <- move$Mo0.1[[3]]}
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
     out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=post[1,], xp0=xp0,
          pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model, Data=Data,
          Status=Status, Thinning=Thinning, Acceptance=Acceptance, Dev=Dev,
          Mon=Mon, Mo0=Mo0, post=post, thinned=thinned, VarCov=VarCov)
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
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
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
