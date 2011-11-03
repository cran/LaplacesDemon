###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The purpose of the LaplaceApproximation function is to maximize the     #
# logarithm of the unnormalized joint posterior distribution of a         #
# Bayesian model with a deterministic approximation that is a conjugate   #
# gradient or adaptive gradient ascent algorithm to estimate the          #
# posterior modes and variances. The gradient ascent algorithm uses an    #
# adaptive step size.                                                     #
###########################################################################

LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="CG", Stop.Tolerance=1.0E-5)
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     if(missing(Model))
          stop("Model is a required argument.")
     if(missing(Data))
          stop("Data is a required argument.")
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
     if((Method != "AGA") & (Method != "CG"))
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
     parm <- as.vector(parm)
     parm.len <- length(parm)
     ans.len <- length(Model(parm, Data)[[1]])
     if(ans.len > 1) stop("Multiple joint posteriors exist!")
     parm.int <- {diag(parm.len) * Interval} / 2
     parm.int <- ifelse(!is.finite(parm.int), 0, parm.int)
     approx.grad <- array(0, dim=c(ans.len, parm.len))
     parm.old <- parm
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- Data$parm.names
     tol.new <- tol.old <- sqrt(sum({parm.new - parm.old}^2))
     Step.Size  <- 0.001
     post <- matrix(parm.new, 1, parm.len)
     iter <- 0
     ### Attempt prior to starting...
     m.old <- Model(parm, Data)
     if(!is.finite(m.old[[1]])) {
          stop("The posterior is not real or numeric.")}
     if(!is.finite(m.old[[2]])) {
          stop("The deviance is not real or numeric.")}
     Dev <- matrix(m.old[[2]],1,1)
     ### Parms could have been constrained...
     parm.old <- m.old[[5]]
     parm.new <- parm.old - 0.1
     tol.new <- tol.old <- sqrt(sum({parm.new - parm.old}^2))
     post <- matrix(parm.new, 1, parm.len)
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     if(Method == "AGA") {
     while({iter < Iterations} & {tol.new > Stop.Tolerance}) {
          iter <- iter + 1
          tol.old <- tol.new
          parm.old <- parm.new
          ### Print Status
          if(iter %% round(Iterations / 10) == 0) {
               cat("Iteration: ", iter, " of ", Iterations, "\n")}
          ### Approximate Truncated Gradient
          for (i in 1:parm.len)
               {
               parm.old <- Model(parm.old, Data)[[5]]
               high <- Model(parm.old + parm.int[,i], Data)[[1]] * -1
               low <- Model(parm.old - parm.int[,i], Data)[[1]] * -1
               approx.grad[,i] <- {high - low} / Interval
               approx.grad <- ifelse(approx.grad > 1000, 1000,
                    approx.grad)
               approx.grad <- ifelse(approx.grad < -1000, -1000,
                    approx.grad)
               if(!is.finite(high) | !is.finite(low)) {
                    approx.grad[,i] <- 0}
               }
          ### In iteration 1 or at status, try 10 different step sizes
          if({iter == 1} | {iter %% round(Iterations / 10) == 0}) {
               step.size <- c(1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3,
                    5.0E-4, 1.0E-4, 5.0E-5, 1.0E-5, 1.0E-6)
               lujpd <- rep(0,10)
               parm.temp <- matrix(0, 10, parm.len)
               for (j in 1:10) {
                    parm.temp[j,] <- parm.old - step.size[j] * c(approx.grad)
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
          parm.new <- parm.old - Step.Size * c(approx.grad)
          parm.new <- ifelse(!is.finite(parm.new), parm.old, parm.new)
          tol.new <- sqrt(sum({parm.new - parm.old}^2))
          m.old <- Model(parm.old, Data)
          m.new <- Model(parm.new, Data)
          Dev.new <- m.old[[2]]
          if(!is.finite(m.new[[1]])) {
               parm.new <- parm.old
               Dev.new <- m.new[[2]]
               tol.new <- 0}
          parm.new <- m.new[[5]]
          post <- rbind(post, parm.new)
          Dev <- rbind(Dev, Dev.new)
          ### Adaptive Step Size
          if(m.new[[1]] <= m.old[[1]]) {
               Step.Size <- Step.Size * 0.999
               if(Step.Size < 1.0E-8) Step.Size <- 1.0E-8
               parm.new <- parm.old}
          }
     ### Logarithm of the Marginal Likelihood
     LML <- LML(Model, Data, parm.new, method="LME2")
     ### Column names to samples
     if(ncol(post) == length(Data$parm.names)) {
          colnames(post) <- Data$parm.names}
     rownames(post) <- 1:{nrow(post)}
     ### Summary
     cat("\nCreating Summary...\n")
     Summ <- matrix(NA, parm.len, 4, dimnames=list(Data$parm.names,
          c("Mode","SD","LB","UB")))
     Summ[,1] <- parm.new
     Summ[,2] <- sqrt(diag(LML[[2]]))
     Summ[,3] <- parm.new - 2*Summ[,2]
     Summ[,4] <- parm.new + 2*Summ[,2]
     time2 <- proc.time()
     ### Output
     converged <- ifelse(tol.new <= Stop.Tolerance, TRUE, FALSE)
     LA <- list(Call=match.call(),
          Converged = converged,
          Covar = LML[[2]],
          Deviance = as.vector(Dev),
          History = post,
          Initial.Values = parm,
          Iterations = iter,
          LML = LML[[1]],
          LP.Final = as.vector(Model(parm.new, Data)[[1]]),
          LP.Initial = Model(parm, Data)[[1]],
          Minutes = round(as.vector(time2[3] - time1[3]) / 60, 2),
          Step.Size.Final = Step.Size,
          Step.Size.Initial = Step.Size.Initial,
          Stop.Tolerance = Stop.Tolerance,
          Summary = Summ,
          Tolerance = sqrt(sum({parm.new - parm.old}^2)))
     }
     if(Method == "CG") {
     ModelWrapper <- function(parm) {
          return(Model(parm, Data)[[1]])
          }
     est <- optim(par=parm, fn=ModelWrapper, gr=NULL,
          method="CG",
          control=list(maxit=Iterations, fnscale=-1), hessian=TRUE)
     ### Hessian
     Inverse.test <- try(VarCov <- -as.inverse(est$hessian), silent=TRUE)
     if(class(Inverse.test) != "try-error") {
          diag(VarCov) <- ifelse(diag(VarCov) <= 0, .Machine$double.eps,
               diag(VarCov))}
     else {
          cat("\nWARNING: Failure to solve matrix inversion of Approx. Hessian.\n")
          cat("NOTE: Identity matrix is supplied instead.\n")
          VarCov <- diag(parm.len)
          }
     ### Logarithm of the Marginal Likelihood
     LML <- NA
     options(warn=-1)
     LML.test <- try(LML <- parm.len/2 * log(2*pi) + 0.5*log(det(VarCov)) +
          as.vector(Model(est$par, Data)[[1]]), silent=TRUE)
     if((class(LML.test) != "try-error") & is.finite(LML.test[1]))
          LML <- LML.test[1]
     options(warn=0)
     ### Summary
     cat("\nCreating Summary...\n")
     Summ <- matrix(NA, parm.len, 4, dimnames=list(Data$parm.names,
          c("Mode","SD","LB","UB")))
     Summ[,1] <- est$par
     Summ[,2] <- sqrt(diag(VarCov))
     Summ[,3] <- est$par - 2*Summ[,2]
     Summ[,4] <- est$par + 2*Summ[,2]
     time2 <- proc.time()
     LA <- list(Call=match.call(),
          Converged=(est$convergence == 0),
          Covar=VarCov,
          Deviance=as.vector(Model(est$par, Data)[[2]]),
          History=NA,
          Initial.Values=parm,
          Iterations=as.vector(est$counts[1]),
          LML=LML[[1]],
          LP.Final=as.vector(Model(est$par, Data)[[1]]),
          LP.Initial=Model(parm, Data)[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Step.Size.Final=NA,
          Step.Size.Initial=NA,
          Stop.Tolerance=NA,
          Summary=Summ,
          Tolerance=NA)
     }
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
     }

#End
