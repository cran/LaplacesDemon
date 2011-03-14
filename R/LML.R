###########################################################################
# LML (Logarithm of the Marginal Likelihood)                              #
#                                                                         #
# This function is hidden from users because it is used only by two       #
# functions: LaplaceApproximation and LaplacesDemon.                      #
###########################################################################

LML <- function(Model, Data, Modes)
     {
     Interval <- 1.0E-6
     parm.len <- length(Modes)
     ### Approximate Hessian
     eps <- Interval * Modes
     Approx.Hessian <- matrix(0, parm.len, parm.len)
     for (i in 1:parm.len) {
          for (j in 1:parm.len) {
               x1 <- Modes
               x1[i] <- x1[i] + eps[i]
               x1[j] <- x1[j] + eps[j]
               x2 <- Modes
               x2[i] <- x2[i] + eps[i]
               x2[j] <- x2[j] - eps[j]
               x3 <- Modes
               x3[i] <- x3[i] - eps[i]
               x3[j] <- x3[j] + eps[j]
               x4 <- Modes
               x4[i] <- x4[i] - eps[i]
               x4[j] <- x4[j] - eps[j]
               Approx.Hessian[i, j] <- (Model(x1, Data)[[1]] -
                    Model(x2, Data)[[1]] - Model(x3, Data)[[1]] +
                    Model(x4, Data)[[1]]) / (4 * eps[i] * eps[j])
               }
          }
     Inverse.test <- try(VarCov.t <- -solve(Approx.Hessian), silent=TRUE)
     if(is.numeric(Inverse.test[1])) {
          VarCov <- -solve(Approx.Hessian)
          diag(VarCov) <- ifelse(diag(VarCov) < 0, diag(VarCov) * -1,
               diag(VarCov)) #Better method of preventing negative variance?
          diag(VarCov) <- ifelse(diag(VarCov) == 0, 1.0E-10, diag(VarCov))}
     if(is.character(Inverse.test[1])) {
          cat("\nWARNING: Failure to solve matrix inversion of Approx. Hessian.\n")
          cat("NOTE: Identity matrix is supplied instead.\n")
          VarCov <- matrix(0, length(Modes), length(Modes))
          diag(VarCov) <- 1
          }
     ### Logarithm of the Marginal Likelihood
     LML <- NA
     options(warn=-1)
     LML.test <- try(LML <- parm.len/2 * log(2*pi) + 0.5*log(det(VarCov)) +
          as.vector(Model(Modes, Data)[[1]]), silent=TRUE)
     if(is.numeric(LML.test[1]) & !is.nan(LML.test[1]) &
          !is.infinite(LML.test[1])) {LML <- LML.test[1]}
     options(warn=0)
     ### Output
     LML.out <- list(LML=LML, VarCov=VarCov)
     return(LML.out)
     }

#End
