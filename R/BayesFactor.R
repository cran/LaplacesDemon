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
     if(class(x[[1]]) != "laplace")
          stop("x is not of class laplace in BayesFactor().")
     if(class(x[[2]]) != "laplace")
          stop("x is not of class laplace in BayesFactor().")
     if(x[[1]]$Converged == FALSE) 
          cat("WARNING: LaplaceApproximation() did not converge in M[1].\n")
     if(x[[2]]$Converged == FALSE) 
          cat("WARNING: LaplaceApproximation() did not converge in M[2].\n")
     LML1 <- x[[1]]$LML
     LML2 <- x[[2]]$LML
     if(is.na(LML1)) stop("LML is missing for M[1].")
     if(is.na(LML2)) stop("LML is missing for M[2].")
     ### Bayes factor
     B <- exp(LML1 - LML2)
     ### Output
     hypothesis <- "M[1] > M[2]"
     if(B <= 0.1) strength <- "Strong against M[1]"
     if((B > 0.1) & (B <= (1/3))) strength <- "Substantial against M[1]"
     if((B > (1/3)) & (B < 1)) {
          strength <- "Barely worth mentioning against M[1]"}
     if((B >= 1) & (B < 3)) strength <- "Barely worth mentioning for M[1]"
     if((B >= 3) & (B < 10)) strength <- "Substantial for M[1]"
     if(B >= 10) strength <- "Strong for M[1]"
     BF.out <- list(B=B, Hypothesis=hypothesis, 
          Strength.of.Evidence=strength)
     return(BF.out)
     }

#End
