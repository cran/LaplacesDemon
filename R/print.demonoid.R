###########################################################################
# print.demonoid                                                          #
#                                                                         #
# The purpose of this function is to print the contents of an object of   #
# class demonoid to the screen.                                           #
###########################################################################

print.demonoid <- function(x, ...)
     {
     if(is.null(x)) stop("x is NULL.\n")
     cat("Call:\n")
     print(x$Call)
     cat("\nAcceptance Rate: ", round(x$Acceptance.Rate,3),
          "\n", sep="")
     cat("Adaptive: ", x$Adaptive, "\n", sep="")
     cat("Algorithm: ", x$Algorithm, "\n", sep="")
     cat("Covar: (NOT SHOWN HERE)\n")
     cat("DIC of all samples (Dbar): ", round(x$DIC1[1],3),
          "\n", sep="")
     cat("DIC of all samples (pD): ", round(x$DIC1[2],3),
          "\n", sep="")
     cat("DIC of all samples (DIC): ", round(x$DIC1[3],3),
          "\n", sep="")
     cat("DIC of stationary samples (Dbar): ",
          round(x$DIC2[1],3), "\n", sep="")
     cat("DIC of stationary samples (pD): ",
          round(x$DIC2[2],3), "\n", sep="")
     cat("DIC of stationary samples (DIC): ",
          round(x$DIC2[3],3), "\n", sep="")
     cat("DR: ", x$DR, "\n", sep="")
     cat("Initial Values:\n")
     print(x$Initial.Values)
     cat("\nIterations: ", x$Iterations, "\n", sep="")
     cat("Minutes of run-time: ", round(x$Minutes,2), "\n",
          sep="")
     cat("Model: (NOT SHOWN HERE)\n")
     cat("Monitor: (NOT SHOWN HERE)\n")
     cat("Parameters (Number of): ", x$Parameters, "\n",
          sep="")
     cat("Periodicity: ", x$Periodicity, "\n", sep="")
     cat("Posterior1: (NOT SHOWN HERE)\n")
     cat("Posterior2: (NOT SHOWN HERE)\n")
     cat("Recommended Burn-In of Thinned Samples: ",
          x$Rec.BurnIn.Thinned, "\n", sep="")
     cat("Recommended Burn-In of Un-thinned Samples: ",
          x$Rec.BurnIn.UnThinned, "\n", sep="")
     cat("Recommended Thinning: ", x$Rec.Thinning, "\n", sep="")
     cat("Status is displayed every ", x$Status,
          " iterations\n", sep="")
     cat("Summary1: (SHOWN BELOW)\n")
     cat("Summary2: (SHOWN BELOW)\n")
     cat("Thinned Samples: ", x$Thinned.Samples, "\n",
          sep="")
     cat("Thinning: ", x$Thinning, "\n", sep="")
     cat("\n\nSummary of All Samples\n")
     print(x$Summary1)
     cat("\n\nSummary of Stationary Samples\n")
     print(x$Summary2)
     invisible(x)
     }

#End
