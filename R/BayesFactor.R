###########################################################################
# BayesFactor                                                             #
#                                                                         #
# The purpose of this function is to estimate a Bayes factor from two     #
# objects of class laplace.                                               #
###########################################################################

BayesFactor <- function(x)
     {
     ### Initial Checks
     if(is.null(x)) stop("x is required in BayesFactor().")
     Model.num <- length(x)
     for (i in 1:Model.num) {
          if((class(x[[i]]) != "demonoid") & (class(x[[i]]) != "laplace"))
            stop("x is not of class demonoid or laplace in BayesFactor().")
          if((class(x[[i]]) == "laplace") && (x[[i]]$Converged == FALSE)) { 
               cat("WARNING: LaplaceApproximation() did not converge in ",
                    "M[",i,"].\n", sep="")}
          if(is.na(x[[i]]$LML))
               stop(cat("LML is missing in M[",i,"].", sep=""))
          }
     ### Bayes factor
     B <- matrix(NA, Model.num, Model.num)
     for (i in 1:Model.num) {for (j in 1:Model.num) {
          B[i,j] <- exp(x[[i]]$LML - x[[j]]$LML)
          }}
     strength <- rep(NA,6)
     strength[1] <- "-Inf  <  B <= 0.1   Strong against"
     strength[2] <- "0.1   <  B <= (1/3) Substantial against"
     strength[3] <- "(1/3) <  B < 1      Barely worth mentioning against"
     strength[4] <- "1     <= B < 3      Barely worth mentioning for"
     strength[5] <- "3     <= B < 10     Substantial for"
     strength[6] <- "10    <= B < Inf    Strong for"
     ### Output
     BF.out <- list(B=B, Hypothesis="row > column", 
          Strength.of.Evidence=strength)
     return(BF.out)
     }

#End
