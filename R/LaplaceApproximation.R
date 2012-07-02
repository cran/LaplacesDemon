###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The purpose of the LaplaceApproximation function is to maximize the     #
# logarithm of the unnormalized joint posterior distribution of a         #
# Bayesian model with a deterministic approximation that is a resilient   #
# backpropagation algorithm (among other selectable algorithms) to        #
# estimate the posterior modes and variances.                             #
###########################################################################

LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="LBFGS", Samples=1000, sir=TRUE,
     Stop.Tolerance=1.0E-5)
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     LA.call <- match.call()
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data$parm.names))}
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n")}}}}
     if({Interval <= 0} | {Interval > 1}) Interval <- 1.0E-6
     Iterations <- round(Iterations)
     if(Iterations < 10) {Iterations <- 10}
     else if(Iterations > 1000000) {Iterations <- 1000000}
     if({Method != "AGA"} & {Method != "CG"} & {Method != "LBFGS"} &
          {Method != "Rprop"})
          stop("Method is unknown.")
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out,1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of apply functions\n")
          cat(     "were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are 'vectorized'.\n")
     }
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are 'vectorized'.\n")
     }
     ### Sample Size of Data
     if(!is.null(Data$n)) if(length(Data$n) == 1) N <- Data$n
     if(!is.null(Data$N)) if(length(Data$N) == 1) N <- Data$N
     if(!is.null(Data$y)) N <- nrow(matrix(Data$y))
     if(!is.null(Data$Y)) N <- nrow(matrix(Data$Y))
     if(!is.null(N)) cat("Sample Size: ", N, "\n")
     else stop("Sample size of Data not found in n, N, y, or Y.")
     ###########################  Preparation  ############################
     ### Attempt prior to starting...
     m.old <- Model(parm, Data)
     if(!is.list(m.old)) stop("Model must return a list.")
     if(length(m.old) != 5) stop("Model must return five components.")
     if(length(m.old[[1]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(parm), length(m.old[[5]])))
          stop("The number of initial values and parameters differs.")
     if(!is.finite(m.old[[1]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          Initial.Values <- GIV(Model, Data)
          m.old <- Model(Initial.Values, Data)
          }
     if(!is.finite(m.old[[1]]))
          stop("The posterior is non-finite.")
     if(!is.finite(m.old[[2]]))
          stop("The deviance is non-finite.")
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     if(Method == "AGA") {
          LA <- AGA(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, LA.call, m.old)}
     if(Method == "CG") {
          cat("WARNING: CG has been deprecated. Rprop is used instead.")
          Method <- "Rprop"}
     if(Method == "LBFGS") {
          LA <- LBFGS(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, LA.call, m.old)}
     if(Method == "Rprop") {
          LA <- Rprop(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, LA.call, m.old)}
     Dev <- as.vector(LA$Dev)
     iter <- LA$iter
     parm.len <- LA$parm.len
     parm.new <- LA$parm.new
     parm.old <- LA$parm.old
     post <- LA$post
     Step.Size <- LA$Step.Size
     tol.new <- LA$tol.new
     if(iter == 1) stop("LaplaceApproximation stopped at iteration 1.")
     converged <- ifelse(tol.new <= Stop.Tolerance, TRUE, FALSE)
     ### Column names to samples
     if(ncol(post) == length(Data$parm.names))
          colnames(post) <- Data$parm.names
     rownames(post) <- 1:{nrow(post)}
     ########################  Covariance Matirx  #########################
     Approx.Hessian <- Hessian(Model, parm.new, Data)
     Inverse.test <- try(-as.inverse(Approx.Hessian), silent=TRUE)
     if(!inherits(Inverse.test, "try-error")) {
          VarCov <- Inverse.test
          diag(VarCov) <- ifelse(diag(VarCov) <= 0,
               .Machine$double.eps, diag(VarCov))}
     else {
          cat("\nWARNING: Failure to solve matrix inversion of ",
               "Approx. Hessian.\n", sep="")
          cat("NOTE: Identity matrix is supplied instead.\n")
               VarCov <- diag(parm.len)}
     #################  Sampling Importance Resampling  ##################
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Sampling from Posterior with Sampling Importance Resampling\n")
          posterior <- SIR(Model, Data, mu=parm.new, Sigma=VarCov,
               n=Samples)
          Mon <- matrix(0, nrow(posterior), length(Data$mon.names))
          dev <- rep(0, nrow(posterior))
          for (i in 1:nrow(posterior)) {
               mod <- Model(posterior[i,], Data)
               dev[i] <- mod[[2]]
               Mon[i,] <- mod[[3]]
               }
          colnames(Mon) <- Data$mon.names}
     else {
          if({sir == TRUE} & {converged == FALSE})
               cat("Posterior samples are not drawn due to Converge=FALSE\n")
          posterior <- NA; Mon <- NA}
     #####################  Summary, Point-Estimate  ######################
     cat("Creating Summary from Point-Estimates\n")
     Summ1 <- matrix(NA, parm.len, 4, dimnames=list(Data$parm.names,
          c("Mode","SD","LB","UB")))
     Summ1[,1] <- parm.new
     Summ1[,2] <- sqrt(diag(VarCov))
     Summ1[,3] <- parm.new - 2*Summ1[,2]
     Summ1[,4] <- parm.new + 2*Summ1[,2]
     ###################  Summary, Posterior Samples  ####################
     Summ2 <- NA
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Creating Summary from Posterior Samples\n")
          Summ2 <- matrix(NA, ncol(posterior), 7,
               dimnames=list(Data$parm.names,
                    c("Mean","SD","MCSE","ESS","LB","Median","UB")))
          Summ2[,1] <- colMeans(posterior)
          Summ2[,2] <- apply(posterior, 2, sd)
          Summ2[,3] <- Summ2[,2] / sqrt(nrow(posterior))
          Summ2[,4] <- rep(nrow(posterior), ncol(posterior))
          Summ2[,5] <- apply(posterior, 2, quantile, c(0.025))
          Summ2[,6] <- apply(posterior, 2, quantile, c(0.500))
          Summ2[,7] <- apply(posterior, 2, quantile, c(0.975))
          Deviance <- rep(0, 7)
          Deviance[1] <- mean(dev)
          Deviance[2] <- sd(dev)
          Deviance[3] <- sd(dev) / sqrt(nrow(posterior))
          Deviance[4] <- nrow(posterior)
          Deviance[5] <- as.numeric(quantile(dev, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(dev, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(dev, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:ncol(Mon)) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(as.vector(Mon[,j]))
               Monitor[3] <- sd(as.vector(Mon[,j])) / sqrt(nrow(Mon))
               Monitor[4] <- nrow(Mon)
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data$mon.names[j]
               }
          }
     ###############  Logarithm of the Marginal Likelihood  ###############
     LML <- LML(Model, Data, parm.new, method="LME2")
     LML <- list(LML=NA, VarCov=VarCov)
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          lml <- LML(theta=posterior, LL=(dev*(-1/2)), method="NSIS")
          LML[[1]] <- lml[[1]]}
     else if({sir == FALSE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(Model, Data, Modes, method="LME2")}
     time2 <- proc.time()
     #############################  Output  ##############################
     LA <- list(Call=LA.call,
          Converged = converged,
          Covar = VarCov,
          Deviance = as.vector(Dev),
          History = post,
          Initial.Values = parm,
          Iterations = iter,
          LML = LML[[1]],
          LP.Final = as.vector(Model(parm.new, Data)[[1]]),
          LP.Initial = Model(parm, Data)[[1]],
          Minutes = round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor = Mon,
          Posterior = posterior,
          Step.Size.Final = Step.Size,
          Step.Size.Initial = 1,
          Summary1 = Summ1,
          Summary2 = Summ2,
          Tolerance.Final = sqrt(sum({parm.new - parm.old}^2)),
          Tolerance.Stop = Stop.Tolerance)
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
     }
AGA <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     LA.call, m.old)
     {
     Dev <- matrix(m.old[[2]],1,1)
     parm <- as.vector(parm)
     parm.len <- length(parm)
     parm.old <- m.old[[5]]
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- Data$parm.names
     tol.new <- tol.old <- sqrt(sum({parm.new - parm.old}^2))
     step.size <- c(1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3,
          5.0E-4, 1.0E-4, 5.0E-5, 1.0E-5, 1.0E-6)
     Step.Size  <- 0.001
     post <- matrix(parm.new, 1, parm.len)
     for (iter in 1:Iterations) {
          tol.old <- tol.new
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Approximate Truncated Gradient
          approx.grad <- partial(Model, parm.old, Data, Interval)
          approx.grad <- interval(approx.grad, -1000, 1000)
          ### In iteration 1 or at status, try 10 different step sizes
          if({iter == 1} | {iter %% round(Iterations / 10) == 0}) {
               lujpd <- rep(0,10)
               parm.temp <- matrix(0, 10, parm.len)
               for (j in 1:10) {
                    parm.temp[j,] <- parm.old + step.size[j] * approx.grad
                    parm.temp[j,] <- ifelse(!is.finite(parm.temp[j,]),
                         parm.old, parm.temp[j,])
                    lujpd[j] <- Model(parm.temp[j,], Data)[[1]]
                    }
               lujpd <- ifelse(!is.finite(lujpd), -1.0E+200, lujpd)
               if(max(lujpd) > m.old[[1]]) {
                    Step.Size <- step.size[which.max(lujpd)]
                    parm.new <- parm.temp[which.max(lujpd),]}
               else {parm.new <- parm.old}
               if(iter == 1) Step.Size.Initial <- Step.Size
               }
          ### Otherwise, work with the current and adaptive Step.Size
          parm.new <- parm.old + Step.Size * approx.grad
          parm.new <- ifelse(!is.finite(parm.new), parm.old, parm.new)
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({m.new[[5]] - parm.old}^2))
          if(!is.finite(m.new[[1]])) {
               parm.new <- parm.old
               tol.new <- 0}
          parm.new <- m.new[[5]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[[2]])
          ### Adaptive Step Size
          if(m.new[[1]] <= m.old[[1]]) {
               Step.Size <- Step.Size * 0.999
               if(Step.Size < 1.0E-8) Step.Size <- 1.0E-8
               parm.new <- parm.old
               m.new <- m.old}
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
LBFGS <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     LA.call, m.old)
     {
     Dev <- matrix(m.old[[2]],1,1)
     parm <- as.vector(parm)
     parm.len <- length(parm)
     parm.new <- parm.old <- m.old[[5]]
     names(parm.new) <- Data$parm.names
     tol <- sqrt(sum({parm.new - parm.old}^2))
     post <- matrix(parm.new, 1, parm.len)
     ModelWrapper <- function(parm.new) {
          out <- Model(parm.new, Data)[[1]]
          return(out)
          }
     for (iter in 1:Iterations)
          {
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### LBFGS
          Fit <- optim(par=parm.new, fn=ModelWrapper,
               method="L-BFGS-B", control=list(fnscale=-1, maxit=1))
          m.new <- Model(Fit$par, Data)
          tol <- sqrt(sum({m.new[[5]] - parm.old}^2))
          if(!is.finite(m.new[[1]]) | {m.new[[1]] < m.old[[1]]}) {
               m.new <- m.old}
          parm.new <- m.new[[5]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[[2]])
          if(tol <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=0,
          tol.new=tol)
     return(LA)
     }
Rprop <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     LA.call, m.old)
     {
     Dev <- matrix(m.old[[2]],1,1)
     parm <- as.vector(parm)
     parm.len <- length(parm)
     approx.grad.old <- approx.grad.new <- rep(0, parm.len)
     parm.old <- m.old[[5]]
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- Data$parm.names
     tol.new <- tol.old <- sqrt(sum({parm.new - parm.old}^2))
     post <- matrix(parm.new, 1, parm.len)
     Step.Size <- rep(0.0125, parm.len)
     for (iter in 1:Iterations)
          {
          approx.grad.old <- approx.grad.new
          tol.old <- tol.new
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Approximate Truncated Gradient
          approx.grad.new <- partial(Model, parm.old, Data, Interval)
          approx.grad.new <- interval(approx.grad.new, -1000, 1000)
          ### Determine if Gradients Changed Sign
          change <- (approx.grad.old >= 0) == (approx.grad.new >= 0)
          ### Adjust weight (step size) based on sign change
          Step.Size[change == TRUE] <- Step.Size[change == TRUE] * 0.5
          Step.Size[change == FALSE] <- Step.Size[change == FALSE] * 1.2
          Step.Size <- interval(Step.Size, 0.0001, 50)
          ### Propose new state based on weighted approximate gradient
          parm.new <- parm.old + (Step.Size * approx.grad.new)
          m.new <- Model(parm.new, Data)
          tol.new <- sqrt(sum({m.new[[5]] - parm.old}^2))
          if(!is.finite(m.new[[1]]) | {m.new[[1]] < m.old[[1]]}) {
               ### New code for testing
               p.order <- sample(1:length(parm.new))
               parm.temp <- parm.old
               for (i in 1:length(p.order)) {
                    parm.temp[p.order[i]] <- parm.new[p.order[i]]
                    m.new <- Model(parm.temp, Data)
                    if(m.new[[1]] < m.old[[1]])
                         parm.temp[p.order[i]] <- parm.old[p.order[i]]}
               m.new <- Model(parm.temp, Data)}
          parm.new <- m.new[[5]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, m.new[[2]])
          if(tol.new <= Stop.Tolerance) break
          }
     Dev <- Dev[-1,]; post <- post[-1,]
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
          parm.old=parm.old, post=post, Step.Size=Step.Size,
          tol.new=tol.new)
     return(LA)
     }
#End
