###########################################################################
# IVO                                                                     #
#                                                                         #
# The Initial Value Optimization (IVO) function maximizes the logarithm   #
# of the unnormalized joint posterior distribution of a Bayesian model    #
# with a deterministic approximation gradient ascent algorithm. This      #
# gradient ascent algorithm uses an adaptive step size.                   #
###########################################################################

IVO <- function(Model=NULL, parm=NULL, Data=NULL, Interval=1.0E-6,
     Iterations=100, Stop.Tolerance=1.0E-5)
     {
     ### Initial Checks
     if(is.null(Model)) stop("Model is a required argument in IVO.\n")
     if(is.null(Data)) stop("Data is a required argument in IVO.\n")
     if(is.null(parm)) {
          cat("Initial values have been set to zero prior to IVO.\n")
          parm <- rep(0, length(Data$Parameters))}
     if((Interval <= 0) | (Interval > 1)) Interval <- 1.0E-6
     Iterations <- round(Iterations)
     if(Iterations < 1) Iterations <- 1
     if(Iterations > 1000000) Iterations <- 1000000
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
     ### Preparation
     parm <- as.vector(parm)
     parm.len <- length(parm)
     ans.len <- length(Model(parm, Data)[[1]])
     parm.int <- (diag(parm.len) * Interval) / 2
     parm.int <- ifelse(is.na(parm.int), 0, parm.int)
     parm.int <- ifelse(is.nan(parm.int), 0, parm.int)
     parm.int <- ifelse(is.infinite(parm.int), 0, parm.int)
     approx.grad <- array(0, dim=c(ans.len, parm.len))
     parm.old <- parm
     parm.new <- parm.old - 0.1 #First step
     tol.new <- tol.old <- sqrt(sum((parm.new - parm.old)^2))
     Step.Size.Initial <- Step.Size  <- 0.001
     iter <- 0
     ### Attempt prior to starting...
     lujpd.old <- Model(parm, Data)[[1]]
     if(is.na(lujpd.old) | is.nan(lujpd.old) | is.infinite(lujpd.old)) {
          stop("The posterior is not real or numeric.")}
     ### Adaptive Gradient Ascent Algorithm
     cat("Initial Values Optimization (IVO) begins...\n")
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
          ### In iteration 1 or at status, try a few different step sizes
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
               if(max(lujpd) > lujpd.old) {
                    Step.Size <- step.size[which.max(lujpd)]
                    parm.new <- parm.temp[which.max(lujpd),]}
               if(max(lujpd) <= lujpd.old) {parm.new <- parm.old}
               }
          ### Otherwise, work with the current and adaptive Step.Size
          parm.new <- parm.old - Step.Size * c(approx.grad)
          parm.new <- ifelse(is.na(parm.new), parm.old, parm.new)
          parm.new <- ifelse(is.nan(parm.new), parm.old, parm.new)
          parm.new <- ifelse(is.infinite(parm.new), parm.old, parm.new)
          tol.new <- sqrt(sum((parm.new - parm.old)^2))
          lujpd.old <- Model(parm.old, Data)[[1]]
          lujpd.new <- Model(parm.new, Data)[[1]]
          if(is.na(lujpd.new) | is.nan(lujpd.new) |
               is.infinite(lujpd.new)) {
               parm.new <- parm.old
               tol.new <- 0}
          ### Adaptive Step Size
          if(lujpd.new <= lujpd.old) {
               Step.Size <- Step.Size * 0.999
               if(Step.Size < 1.0E-8) Step.Size <- 1.0E-8
               parm.new <- parm.old}
          }
     ### Output
     converged <- ifelse(tol.new <= Stop.Tolerance, TRUE, FALSE)
     IVO <- list(Converged = converged,
          Iterations=iter,
          LP.Final = Model(parm.new, Data)[[1]],
          LP.Initial = Model(parm, Data)[[1]],
          Parms.Final = parm.new,
          Parms.Initial = parm,
          Step.Size.Final = Step.Size,
          Step.Size.Initial = Step.Size.Initial,
          Stop.Tolerance = Stop.Tolerance,
          Tolerance = sqrt(sum((parm.new - parm.old)^2)))
     cat("Initial Values Optimization (IVO) is finished.\n")
     return(IVO)
     }

#End
