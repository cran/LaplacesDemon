###########################################################################
# print.laplace                                                           #
#                                                                         #
# The purpose of this function is to print the contents of an object of   #
# class laplace to the screen.                                            #
###########################################################################

print.laplace <- function(x, ...)
     {
     if(is.null(x)) stop("x is NULL.\n")
     cat("\nCall:\n")
     print(x$Call)
     cat("\nConverged: ", x$Converged, "\n", sep="")
     cat("Covariance Matrix: (NOT SHOWN HERE; diagonal shown instead)\n")
     print(diag(x$Covar))
     cat("\nDeviance (Final): ", x$Deviance[length(x$Deviance)], "\n")
     cat("History: (NOT SHOWN HERE)\n")
     cat("Initial Values:\n")
     print(x$Initial.Values)
     cat("\nIterations: ", x$Iterations, "\n", sep="")
     cat("Log(Marginal Likelihood): ", x$LML, "\n", sep="")
     cat("Log-Posterior (Final): ", x$LP.Final, "\n", sep="")
     cat("Log-Posterior (Initial): ", x$LP.Initial, "\n", sep="")
     cat("Minutes of run-time: ", x$Minutes, "\n", sep="")
     cat("Step Size (Final): ", x$Step.Size.Final, "\n", sep="")
     cat("Step Size (Initial): ", x$Step.Size.Initial, "\n", sep="")
     cat("Stop Tolerance: ", x$Stop.Tolerance, "\n", sep="")
     cat("Summary: (SHOWN BELOW)\n")
     cat("Tolerance: ", x$Tolerance, "\n", sep="")
     cat("\nSummary:\n")
     print(x$Summary)
     invisible(x)
     }

#End
