###########################################################################
# Consort                                                                 #
#                                                                         #
# The purpose of this function is to consort with Laplace's Demon         #
# regarding an object of class demonoid.                                  #
###########################################################################

Consort <- function(object=NULL)
     {
     if(is.null(object)) {cat("ERROR: The object argument is empty.\n")}
     oname <- deparse(substitute(object))
     dname <- as.vector(strsplit(as.character(object$Call), "=")[3][[1]])
     cat("\n#############################################################\n")
     cat("# Consort with Laplace's Demon                              #\n")
     cat("#############################################################\n")
     print.demonoid(object)
     ### Check Acceptance.Rate
     Acc.Rate.Level <- 2
     Acc.Rate.Low <- 0.15
     Acc.Rate.High <- 0.5
     if(object$Acceptance.Rate < Acc.Rate.Low) Acc.Rate.Level <- 1
     if(object$Acceptance.Rate > Acc.Rate.High) Acc.Rate.Level <- 3
     ### Check MCSE
     MCSE.crit <- 0.0627
     MCSE.temp <- object$Summary2[1:object$Parameters,3] /
          object$Summary2[1:object$Parameters,2]
     MCSE.temp2 <- sum(is.na(MCSE.temp))
     MCSE.tot <- 0
     if(MCSE.temp2 < object$Parameters) {
          MCSE.tot <- sum(MCSE.temp < MCSE.crit)}
     ### Check Eff.Size
     Eff.Size.temp <- object$Summary2[1:object$Parameters,4]
     Eff.Size.temp <- ifelse(is.na(Eff.Size.temp), 0, Eff.Size.temp)
     Eff.Size.temp <- ifelse(is.nan(Eff.Size.temp), 0, Eff.Size.temp)
     Eff.Size.temp <- ifelse(is.infinite(Eff.Size.temp), 0, Eff.Size.temp)
     Eff.Size.min <- min(Eff.Size.temp)
     Eff.Size.crit <- 100
     ### Check Stationarity
     Stationarity <- FALSE
     if(object$Rec.BurnIn.Thinned < object$Thinned.Samples) {
          Stationarity <- TRUE}
     ### Suggested Values
     Rec.Iterations <- trunc(object$Rec.Thinning / object$Thinning *
          object$Iterations)
     Rec.Periodicity <- object$Rec.Thinning
     Status.temp <- round(Rec.Iterations / (object$Minutes *
          Rec.Iterations / object$Iterations),0)
     if(Status.temp < Rec.Iterations) Rec.Status <- Status.temp
     else Rec.Status <- sqrt(Rec.Iterations)
     ### The Demonic Suggestion of Laplace's Demon
     cat("\nDemonic Suggestion\n\n")
     
     cat("Due to the combination of the following conditions,\n\n")

     cat("1. ", object$Algorithm, "\n", sep="")
     if(Acc.Rate.Level == 1) {
          cat("2. The acceptance rate (", object$Acceptance.Rate,
               ") is below ", Acc.Rate.Low, ".\n", sep="")}
     if(Acc.Rate.Level == 2) {
          cat("2. The acceptance rate (", object$Acceptance.Rate,
               ") is within the interval of [", Acc.Rate.Low, ",",
               Acc.Rate.High, "].\n", sep="")}
     if(Acc.Rate.Level == 3) {
          cat("2. The acceptance rate (", object$Acceptance.Rate,
               ") is above ", Acc.Rate.High, ".\n", sep="")}
     if(MCSE.tot < object$Parameters) {
          cat("3. At least one target MCSE is >= ",
               MCSE.crit * 100,"% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     if(MCSE.tot >= object$Parameters) {
          cat("3. Each target MCSE is < ", MCSE.crit * 100,
               "% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     if(Eff.Size.min < Eff.Size.crit) {
          cat("4. At least one target distribution has an ",
               "effective sample size\n", sep="")
          cat("   less than ", Eff.Size.crit, ".\n", sep="")}
     if(Eff.Size.min >= Eff.Size.crit) {
          cat("4. Each target distribution has an effective ",
               "sample size of\n", sep="")
          cat("   at least ", Eff.Size.crit, ".\n", sep="")}
     if(Stationarity == FALSE) {
          cat("5. At least one target distribution is not ",
               "stationary.\n\n", sep="")}
     if(Stationarity == TRUE) {
          cat("5. Each target distribution became stationary by\n")
          cat("   ", object$Rec.BurnIn.Thinned,
               " iterations.\n\n", sep="")}

     ### Determine if Laplace's Demon is appeased...
     Appeased <- FALSE
     if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
          Appeased <- TRUE
          cat("Laplace's Demon has been appeased, and suggests\n")
          cat("the marginal posterior samples should be plotted\n")
          cat("and subjected to any other MCMC diagnostic deemed\n")
          cat("fit before using these samples for inference.\n\n")
          }
     if(Appeased == FALSE) {
          cat("Laplace's Demon has not been appeased, and suggests\n")
          cat("copy/pasting the following R code into the R",
               " console,\n", sep="")
          cat("and running it.\n\n")

          cat("Initial.Values <- ", oname, "$Posterior1",
                    "[", oname, "$Thinned.Samples,]\n", sep="")
          ### If adaptive, low acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Log.Posterior, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                   "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, low acc., good MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, good acc., good MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If adaptive, high acc., good MCSE, good Eff, good Stat...
          if((object$Adaptive <= object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ## If nonadapt, low acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, low acc., good MCSE, good Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 1) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### DRM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=1, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat("object <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, good acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 2) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", 0,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., bad MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., good MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., bad MCSE, bad Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., good MCSE, bad Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min < Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., bad MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., good MCSE, good Eff, bad Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == FALSE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., bad MCSE, good Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot < object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          ### If nonadapt, high acc., good MCSE, good Eff, good Stat...
          if((object$Adaptive > object$Iterations) &
          (Acc.Rate.Level == 3) & (MCSE.tot >= object$Parameters) &
          (Eff.Size.min >= Eff.Size.crit) & (Stationarity == TRUE)) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    object$Rec.Thinning, ")\n\n", sep="")
               }
          }
     cat("Laplace's Demon is finished consorting.\n")
     }

#End
