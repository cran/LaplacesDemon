###########################################################################
# Consort                                                                 #
#                                                                         #
# The purpose of the Consort function is to consort with Laplace's Demon  #
# regarding an object of class demonoid.                                  #
###########################################################################

Consort <- function(object=NULL)
     {
     if(is.null(object)) stop("The object argument is empty.")
     if(class(object) != "demonoid")
          stop("Consort requires an object of class demonoid.")
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
     if(any(object$Acceptance.Rate == 0)) {
       cat("\nWARNING: Acceptance Rate = 0\n\n")
       cat(oname, " <- LaplacesDemon(Model, Data=", dname,
            ", Adaptive=0,\n", sep="")
       cat("     Covar=NULL, DR=0, Gibbs=FALSE, Initial.Values, Iterations=1000000,\n",
            sep="")
       cat("     Mixture=TRUE, Periodicity=0, Status=1000, Thinning=1000)\n\n", sep="")
       stop("Try the above code before consorting again.")
       }
     if(any(object$Acceptance.Rate < Acc.Rate.Low)) {Acc.Rate.Level <- 1}
     else if(any(object$Acceptance.Rate > Acc.Rate.High)) {Acc.Rate.Level <- 3}
     ### Check MCSE
     MCSE.crit <- 0.0627
     MCSE.temp <- object$Summary2[1:object$Parameters,3] /
          object$Summary2[1:object$Parameters,2]
     MCSE.temp2 <- sum(!is.finite(MCSE.temp))
     MCSE.temp <- ifelse(is.finite(MCSE.temp), MCSE.temp, 0)
     MCSE.tot <- 0
     if(MCSE.temp2 < object$Parameters) {
          MCSE.tot <- sum(MCSE.temp < MCSE.crit)}
     ### Check ESS
     ESS.temp <- object$Summary2[1:object$Parameters,4]
     ESS.temp <- ifelse(!is.finite(ESS.temp), 0, ESS.temp)
     ESS.min <- min(ESS.temp)
     ESS.crit <- 100
     ### Check Stationarity
     Stationarity <- FALSE
     if(object$Rec.BurnIn.Thinned < object$Thinned.Samples) {
          Stationarity <- TRUE}
     ### Suggested Values
     Rec.Iterations <- trunc(object$Rec.Thinning / object$Thinning *
          object$Iterations)
     if(any(object$Acceptance.Rate == 0)) {
          Rec.Adaptive <- round(object$Parameters * 1/1.0E-7)
          Rec.Periodicity <- round(object$Rec.Thinning * 1/1.0E-7)}
     else if(all(object$Acceptance.Rate > 0)) {
          Rec.Adaptive <- round(object$Parameters *
               1/mean(object$Acceptance.Rate))
          Rec.Periodicity <- round(object$Rec.Thinning *
               1/mean(object$Acceptance.Rate))}
     if(Rec.Adaptive > Rec.Iterations) Rec.Adaptive <- Rec.Iterations - 1
     if(Rec.Periodicity > Rec.Iterations) {
          Rec.Periodicity <- Rec.Iterations - 1}
     Status.temp <- trunc(Rec.Iterations / {object$Minutes *
          Rec.Iterations / object$Iterations})
     if(Status.temp < Rec.Iterations) Rec.Status <- Status.temp
     else Rec.Status <- trunc(sqrt(Rec.Iterations))
     Rec.Thinning <- object$Rec.Thinning
     ### The Demonic Suggestion of Laplace's Demon
     cat("\nDemonic Suggestion\n\n")

     cat("Due to the combination of the following conditions,\n\n")

     cat("1. ", object$Algorithm, "\n", sep="")
     if(Acc.Rate.Level == 1) {
          cat("2. The acceptance rate (", min(object$Acceptance.Rate),
               ") is below ", Acc.Rate.Low, ".\n", sep="")}
     else if(Acc.Rate.Level == 2) {
          cat("2. The acceptance rate (", mean(object$Acceptance.Rate),
               ") is within the interval [", Acc.Rate.Low, ",",
               Acc.Rate.High, "].\n", sep="")}
     else {
          cat("2. The acceptance rate (", max(object$Acceptance.Rate),
               ") is above ", Acc.Rate.High, ".\n", sep="")}
     if(MCSE.tot < object$Parameters) {
          cat("3. At least one target MCSE is >= ",
               MCSE.crit * 100,"% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     else {
          cat("3. Each target MCSE is < ", MCSE.crit * 100,
               "% of its marginal posterior\n", sep="")
          cat("   standard deviation.\n")}
     if(ESS.min < ESS.crit) {
          cat("4. At least one target distribution has an ",
               "effective sample size\n", sep="")
          cat("   (ESS) less than ", ESS.crit, ".\n", sep="")}
     else {
          cat("4. Each target distribution has an effective ",
               "sample size (ESS)\n", sep="")
          cat("   of at least ", ESS.crit, ".\n", sep="")}
     if(Stationarity == FALSE) {
          cat("5. At least one target distribution is not ",
               "stationary.\n\n", sep="")}
     else {
          cat("5. Each target distribution became stationary by\n")
          cat("   ", object$Rec.BurnIn.Thinned,
               " iterations.\n\n", sep="")}

     ### Determine if Laplace's Demon is appeased...
     Appeased <- FALSE
     if({object$Adaptive > object$Iterations} &
          {Acc.Rate.Level == 2} & {MCSE.tot >= object$Parameters} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE}) {
          Appeased <- TRUE
          cat("Laplace's Demon has been appeased, and suggests\n")
          cat("the marginal posterior samples should be plotted\n")
          cat("and subjected to any other MCMC diagnostic deemed\n")
          cat("fit before using these samples for inference.\n\n")
          }
     else {
          cat("Laplace's Demon has not been appeased, and suggests\n")
          cat("copy/pasting the following R code into the R",
               " console,\n", sep="")
          cat("and running it.\n\n")

          cat("Initial.Values <- as.initial.values(", oname, ")\n", sep="")
          if({object$Algorithm == "Adaptive Metropolis-within-Gibbs"} |
              {object$Algorithm == "Metropolis-within-Gibbs"}) 
               Time <- object$Iterations / object$Minutes
          else Time <- object$Iterations / object$Minutes / object$Parameters
          ### If >= 100 componetwise iterations per minute and ready...
          if({Time >= 100} &
          ({Acc.Rate.Level == 2} & {MCSE.tot >= object$Parameters} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE})) {
               ### MWG
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, Gibbs=TRUE, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Mixture=FALSE, Periodicity=0", 
                    ", Status=", Rec.Status, ", Thinning=",
                    Rec.Thinning, ")\n\n", sep="")
               }
          ### If >= 100 componetwise iterations per minute not ready...
          if({Time >= 100} &
          ({Acc.Rate.Level != 2} | {MCSE.tot < object$Parameters} |
          {ESS.min < ESS.crit} | {Stationarity == FALSE})) {
               ### AMWG
               if({object$Algorithm != "Adaptive Metropolis-within-Gibbs"} &
                    {object$Algorithm != "Metropolis-within-Gibbs"}) {
                    Rec.Status <- trunc(Rec.Status / object$Parameters)}
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=2,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, Gibbs=TRUE, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Mixture=FALSE, Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    Rec.Thinning, ")\n\n", sep="")
               }
          ### If < 100 componetwise iterations per minute and ready...
          if({Time < 100} &
          ({Acc.Rate.Level == 2} & {MCSE.tot >= object$Parameters} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE})) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=0,\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, Gibbs=FALSE, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Mixture=TRUE, Periodicity=0", 
                    ", Status=", Rec.Status, ", Thinning=",
                    Rec.Thinning, ")\n\n", sep="")
               }
          ### If < 100 componetwise iterations per minute and not ready...
          if({Time < 100} &
          ({Acc.Rate.Level != 2} | {MCSE.tot < object$Parameters} |
          {ESS.min < ESS.crit} | {Stationarity == FALSE})) {
               ### AMM
               if({object$Algorithm == "Adaptive Metropolis-within-Gibbs"} |
                    {object$Algorithm == "Metropolis-within-Gibbs"}) {
                    Rec.Status <- Rec.Status * object$Parameters}
               cat(oname, " <- LaplacesDemon(Model, ",
                    "Data=", dname, ", Adaptive=", Rec.Adaptive, ",\n", sep="")
               cat("     Covar=", oname, "$Covar, DR=0, Gibbs=FALSE, ",
                    "Initial.Values, Iterations=", Rec.Iterations,
                    ",\n", sep="")
               cat("     Mixture=TRUE, Periodicity=", Rec.Periodicity,
                    ", Status=", Rec.Status, ", Thinning=",
                    Rec.Thinning, ")\n\n", sep="")
               }
          }
     cat("Laplace's Demon is finished consorting.\n")
     }

#End
