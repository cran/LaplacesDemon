###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The LaplaceApproximation function maximizes the logarithm of the        #
# unnormalized joint posterior distribution of a Bayesian model with a    #
# deterministic approximation gradient ascent algorithm to estimate the   #
# posterior modes and variances. This gradient ascent algorithm uses an   #
# adaptive step size.                                                     #
###########################################################################

LaplaceApproximation <- function(Model=NULL, parm=NULL, Data=NULL,
     Interval=1.0E-6, Iterations=100, Stop.Tolerance=1.0E-5)
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     if(is.null(Model))
          stop("Model is a required argument in LaplaceApproximation().\n")
     if(is.null(Data))
          stop("Data is a required argument in LaplaceApproximation().\n")
     if(is.null(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data$Parameters))}
     if((Interval <= 0) | (Interval > 1)) Interval <- 1.0E-6
     Iterations <- round(Iterations)
     if(Iterations < 10) Iterations <- 10
     if(Iterations > 1000000) Iterations <- 1000000
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
     as.character.function <- function(x, ... )
          {
          fname <- deparse( substitute( x ) )
          f <- match.fun( x ) 
          out <- c( sprintf( '"%s" <- ', fname), capture.output( f ) ) 
          if(grepl( "^[<]", tail(out,1))) out <- head( out, -1)
          out
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
     if(is.null(N)) stop("Sample size of Data not found in n, N, y, or Y.")
     ###########################  Preparation  ############################
     parm <- as.vector(parm)
     parm.len <- length(parm)
     ans.len <- length(Model(parm, Data)[[1]])
     if(ans.len > 1) stop("Multiple joint posteriors exist!\n")
     parm.int <- (diag(parm.len) * Interval) / 2
     parm.int <- ifelse(is.na(parm.int), 0, parm.int)
     parm.int <- ifelse(is.nan(parm.int), 0, parm.int)
     parm.int <- ifelse(is.infinite(parm.int), 0, parm.int)
     approx.grad <- array(0, dim=c(ans.len, parm.len))
     parm.old <- parm
     parm.new <- parm.old - 0.1 #First step
     names(parm.new) <- Data$parm.names
     tol.new <- tol.old <- sqrt(sum((parm.new - parm.old)^2))
     Step.Size  <- 0.001
     post <- matrix(parm.new, 1, parm.len)
     iter <- 0
     ### Attempt prior to starting...
     m.old <- Model(parm, Data)
     if(is.na(m.old[[1]]) | is.nan(m.old[[1]]) | is.infinite(m.old[[1]])) {
          stop("The posterior is not real or numeric.")}
     if(is.na(m.old[[2]]) | is.nan(m.old[[2]]) | is.infinite(m.old[[2]])) {
          stop("The deviance is not real or numeric.")}
     Dev <- matrix(m.old[[2]],1,1)
     ### Parms could have been constrained...
     parm.old <- m.old[[5]]
     parm.new <- parm.old - 0.1
     tol.new <- tol.old <- sqrt(sum((parm.new - parm.old)^2))
     post <- matrix(parm.new, 1, parm.len)
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     while((iter < Iterations) & (tol.new > Stop.Tolerance)) {
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
               approx.grad[,i] <- (high - low) / Interval
               approx.grad <- ifelse(approx.grad > 1000, 1000,
                    approx.grad)
               approx.grad <- ifelse(approx.grad < -1000, -1000,
                    approx.grad)
               if(is.na(high) | is.nan(high) |
                    is.infinite(high) | is.na(low) |
                    is.nan(low) | is.infinite(low)) {
                    approx.grad[,i] <- 0}
               }
          ### In iteration 1 or at status, try 10 different step sizes
          if((iter == 1) | (iter %% round(Iterations / 10) == 0)) {
               step.size <- c(1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3,
                    5.0E-4, 1.0E-4, 5.0E-5, 1.0E-5, 1.0E-6)
               lujpd <- rep(0,10)
               parm.temp <- matrix(0, 10, parm.len)
               for (j in 1:10) {
                    parm.temp[j,] <- parm.old - step.size[j] * c(approx.grad)
                    parm.temp[j,] <- ifelse(is.na(parm.temp[j,]), parm.old,
                         parm.temp[j,])
                    parm.temp[j,] <- ifelse(is.nan(parm.temp[j,]), parm.old,
                         parm.temp[j,])
                    parm.temp[j,] <- ifelse(is.infinite(parm.temp[j,]),
                         parm.old, parm.temp[j,])
                    lujpd[j] <- Model(parm.temp[j,], Data)[[1]]
                    }
               lujpd <- ifelse(is.na(lujpd) | is.nan(lujpd) |
                    is.infinite(lujpd), -1.0E+200, lujpd)
               if(max(lujpd) > m.old[[1]]) {
                    Step.Size <- step.size[which.max(lujpd)]
                    parm.new <- parm.temp[which.max(lujpd),]}
               if(max(lujpd) <= m.old[[1]]) {parm.new <- parm.old}
               if(iter == 1) Step.Size.Initial <- Step.Size
               }
          ### Otherwise, work with the current and adaptive Step.Size
          parm.new <- parm.old - Step.Size * c(approx.grad)
          parm.new <- ifelse(is.na(parm.new), parm.old, parm.new)
          parm.new <- ifelse(is.nan(parm.new), parm.old, parm.new)
          parm.new <- ifelse(is.infinite(parm.new), parm.old, parm.new)
          tol.new <- sqrt(sum((parm.new - parm.old)^2))
          m.old <- Model(parm.old, Data)
          m.new <- Model(parm.new, Data)
          Dev.new <- m.old[[2]]
          if(is.na(m.new[[1]]) | is.nan(m.new[[1]]) |
               is.infinite(m.new[[1]])) {
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
     LML <- LML(Model, Data, parm.new)
     ### Column names to samples
     if(NCOL(post) == length(Data$parm.names)) {
          colnames(post) <- Data$parm.names}
     rownames(post) <- 1:(NROW(post))
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
          Iterations=iter,
          LML=LML[[1]],
          LP.Final = as.vector(Model(parm.new, Data)[[1]]),
          LP.Initial = Model(parm, Data)[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Step.Size.Final = Step.Size,
          Step.Size.Initial = Step.Size.Initial,
          Stop.Tolerance = Stop.Tolerance,
          Summary = Summ,
          Tolerance = sqrt(sum((parm.new - parm.old)^2)))
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
     }

#End
