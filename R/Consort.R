###########################################################################
# Consort                                                                 #
#                                                                         #
# The purpose of the Consort function is to consort with Laplace's Demon  #
# regarding an object of class demonoid.                                  #
###########################################################################

Consort <- function(object=NULL)
     {
     if(is.null(object)) stop("The object argument is empty.")
     if(!identical(class(object), "demonoid"))
          stop("Consort requires an object of class demonoid.")
     oname <- deparse(substitute(object))
     dname <- as.vector(strsplit(as.character(object$Call), "=")[3][[1]])
     cat("\n#############################################################\n")
     cat("# Consort with Laplace's Demon                              #\n")
     cat("#############################################################\n")
     print.demonoid(object)
     ### Check Acceptance.Rate
     Acc.Rate.Level <- 2
     if((object$Algorithm == "Adaptive Hamiltonian Monte Carlo") |
        (object$Algorithm == "Tempered Hamiltonian Monte Carlo")) {
          L <- substr(object$Call, 1, nchar(object$Call))
          L <- strsplit(L, " ")
          L <- as.numeric(sub(",", "",
               L[[length(L)]][[length(L[[length(L)]])-3]]))
          if(L == 1) {
               Acc.Rate.Low <- 0.5
               Acc.Rate.High <- 0.65}
          else {
               Acc.Rate.Low <- 0.6
               Acc.Rate.High <- 0.7}}
     else if((object$Algorithm == "Componentwise Hit-And-Run Metropolis") |
          (object$Algorithm == "Reversible-Jump")) {
          Acc.Rate.Low <- 0.10
          Acc.Rate.High <- 0.90
          }
     else if(object$Algorithm == "Slice Sampler") {
          Acc.Rate.Low <- 1
          Acc.Rate.High <- 1
          }
     else if(object$Algorithm == "Hamiltonian Monte Carlo") {
          L <- substr(object$Call, 1, nchar(object$Call))
          L <- strsplit(L, " ")
          L <- as.numeric(sub(")", "",
                    L[[length(L)]][[length(L[[length(L)]])]]))
          if(L == 1) {
               Acc.Rate.Low <- 0.5
               Acc.Rate.High <- 0.65}
          else {
               Acc.Rate.Low <- 0.6
               Acc.Rate.High <- 0.7}}
     else if(object$Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          A <- substr(object$Call, 1, nchar(object$Call))
          A <- strsplit(A, " ")
          A <- as.numeric(sub(",", "",
               A[[length(A)]][[length(A[[length(A)]])-12]]))
          delta <- substr(object$Call, 1, nchar(object$Call))
          delta <- strsplit(delta, " ")
          delta <- as.numeric(sub(",", "",
               delta[[length(delta)]][[length(delta[[length(delta)]])-9]]))
          Lmax <- substr(object$Call, 1, nchar(object$Call))
          Lmax <- strsplit(Lmax, " ")
          Lmax <- as.numeric(sub(",", "",
               Lmax[[length(Lmax)]][[length(Lmax[[length(Lmax)]])-3]]))
          lambda <- substr(object$Call, 1, nchar(object$Call))
          lambda <- strsplit(lambda, " ")
          lambda <- as.numeric(sub(")", "",
               lambda[[length(lambda)]][[length(lambda[[length(lambda)]])]]))
          Acc.Rate.Low <- max(round(delta - 0.05, 2), 0.01)
          Acc.Rate.High <- min(round(delta + 0.05, 2), 1)}
     else if(object$Algorithm == "No-U-Turn Sampler") {
          A <- substr(object$Call, 1, nchar(object$Call))
          A <- strsplit(A, " ")
          A <- as.numeric(sub(",", "",
               A[[length(A)]][[length(A[[length(A)]])-6]]))
          delta <- substr(object$Call, 1, nchar(object$Call))
          delta <- strsplit(delta, " ")
          delta <- as.numeric(sub(",", "",
               delta[[length(delta)]][[length(delta[[length(delta)]])-3]]))
          Acc.Rate.Low <- max(round(delta - 0.05, 2), 0.01)
          Acc.Rate.High <- min(round(delta + 0.05, 2), 1)}
     else {
          Acc.Rate.Low <- 0.15
          Acc.Rate.High <- 0.5}
     if(any(object$Acceptance.Rate == 0)) {
          cat("\nWARNING: Acceptance Rate = 0\n\n")
          cat(oname, " <- LaplacesDemon(Model, Data=", dname,
               ", Initial.Values,\n", sep="")
          cat("     Covar=NULL, Iterations=100000, Status=1000, ",
               "Thinning=100,\n", sep="")
          cat("     Algorithm=\"AMM\", Specs=list(Adaptive=500, ",
               "Periodicity=100, w=0.05))\n\n", sep="")
          stop("Try the above code before consorting again.")}
     if(any(object$Acceptance.Rate < Acc.Rate.Low)) {Acc.Rate.Level <- 1}
     else if(any(object$Acceptance.Rate > Acc.Rate.High)) {Acc.Rate.Level <- 3}
     ### Check MCSE
     MCSE.crit <- 0.0627
     MCSE.temp <- object$Summary2[1:object$Parameters,3] /
          object$Summary2[1:object$Parameters,2]
     MCSE.temp2 <- sum(!is.finite(MCSE.temp))
     MCSE.temp <- ifelse(is.finite(MCSE.temp), MCSE.temp, 0)
     MCSE.tot <- 0
     if(MCSE.temp2 < object$Parameters)
          MCSE.tot <- sum(MCSE.temp < MCSE.crit)
     ### Check ESS
     ESS.temp <- object$Summary2[1:object$Parameters,4]
     ESS.temp <- ifelse(!is.finite(ESS.temp), 0, ESS.temp)
     ESS.min <- min(ESS.temp)
     ESS.crit <- 100
     ### Check Stationarity
     Stationarity <- FALSE
     if(object$Rec.BurnIn.Thinned < object$Thinned.Samples)
          Stationarity <- TRUE
     ### Check Diminishing Adaptation (If Adaptive)
     if(nrow(object$CovarDHis) > 1) {
          Dim.Adapt <- as.vector(lm(rowMeans(diff(object$CovarDHis)) ~ c(1:nrow(diff(object$CovarDHis))))$coef[2]) <= 0
          }
     else Dim.Adapt <- TRUE
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
     if(Rec.Status > Rec.Iterations) Rec.Status <- Rec.Iterations
     Rec.Thinning <- object$Rec.Thinning
     ### Miscellaneous
     if(object$Algorithm == "Tempered Hamiltonian Monte Carlo") {
          Temperature <- substr(object$Call, 1, nchar(object$Call))
          Temperature <- strsplit(Temperature, " ")
          Temperature <- as.numeric(sub(")", "",
               Temperature[[length(Temperature)]][[length(Temperature[[length(Temperature)]])]]))}
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
          if(Dim.Adapt == FALSE) {
               cat("WARNING: Diminishing adaptation did not occur.\n")
               if({object$Algorithm != "No-U-Turn Sampler"} &
                  {object$Algorithm != "Sequential Adaptive Metropolis-within-Gibbs"} &
                  {object$Algorithm != "Updating Sequential Adaptive Metropolis-within-Gibbs"})
                    cat("         A new algorithm will be suggested.\n\n")}
          
          cat("Laplace's Demon has not been appeased, and suggests\n")
          cat("copy/pasting the following R code into the R",
               " console,\n", sep="")
          cat("and running it.\n\n")

          cat("Initial.Values <- as.initial.values(", oname, ")\n", sep="")
          if({object$Algorithm == "Adaptive Metropolis-within-Gibbs"} |
              {object$Algorithm == "Metropolis-within-Gibbs"}) 
               Time <- object$Iterations / object$Minutes
          else Time <- object$Iterations / object$Minutes / object$Parameters
          if(Time >= 100) Fast <- TRUE
          else Fast <- FALSE

          if({Acc.Rate.Level == 2} & {MCSE.tot >= object$Parameters} &
          {ESS.min >= ESS.crit} & {Stationarity == TRUE}) Ready <- TRUE
          else Ready <- FALSE

          if(object$Algorithm == "Adaptive Hamiltonian Monte Carlo") Alg <- "AHMC"
          else if(object$Algorithm == "Adaptive Metropolis") Alg <- "AM"
          else if(object$Algorithm == "Adaptive-Mixture Metropolis") Alg <- "AMM"
          else if(object$Algorithm == "Adaptive Metropolis-within-Gibbs") Alg <- "AMWG"
          else if(object$Algorithm == "Componentwise Hit-And-Run Metropolis") Alg <- "CHARM"
          else if(object$Algorithm == "Delayed Rejection Adaptive Metropolis") Alg <- "DRAM"
          else if(object$Algorithm == "Delayed Rejection Metropolis") Alg <- "DRM"
          else if(object$Algorithm == "Experimental") Alg <- "Exper"
          else if(object$Algorithm == "Hamiltonian Monte Carlo") Alg <- "HMC"
          else if(object$Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging")
               Alg <- "HMCDA"
          else if(object$Algorithm == "Hit-And-Run Metropolis") Alg <- "HARM"
          else if(object$Algorithm == "Independence Metropolis") Alg <- "IM"
          else if(object$Algorithm == "Interchain Adaptation") Alg <- "INCA"
          else if(object$Algorithm == "Metropolis-within-Gibbs") Alg <- "MWG"
          else if(object$Algorithm == "No-U-Turn Sampler") Alg <- "NUTS"
          else if(object$Algorithm == "Robust Adaptive Metropolis") Alg <- "RAM"
          else if(object$Algorithm == "Random-Walk Metropolis") Alg <- "RWM"
          else if(object$Algorithm == "Reversible-Jump") Alg <- "RJ"
          else if(object$Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") Alg <- "SAMWG"
          else if(object$Algorithm == "Sequential Metropolis-within-Gibbs") Alg <- "SMWG"
          else if(object$Algorithm == "Slice Sampler") Alg <- "Slice"
          else if(object$Algorithm == "Tempered Hamiltonian Monte Carlo") Alg <- "THMC"
          else if(object$Algorithm == "t-walk") Alg <- "t-walk"
          else if(object$Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") Alg <- "USAMWG"
          else Alg <- "USMWG"

          if({(Alg == "AHMC") & Dim.Adapt & !Ready} |
             {(Alg == "HMC") & Dim.Adapt & !Ready}) {
               ### AHMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AHMC\",\n", sep="")
               cat("     Specs=list(epsilon=", oname, "$CovarDHis[nrow(",
                    oname, "$CovarDHis),],\n", sep="")
               cat("     L=", L, ", Periodicity=", object$Periodicity,
                    "))\n\n", sep="")
               }
          if({(Alg == "AM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "AM") & Dim.Adapt & !Fast & !Ready}) {
               ### AM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "))\n\n", sep="")
               }
          if({(Alg == "AHMC") & !Dim.Adapt} |
             {(Alg == "AM") & !Dim.Adapt & !Fast & Ready} |
             {(Alg == "AMM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "AMM") & Dim.Adapt & !Fast & !Ready} |
             {(Alg == "DRAM") & !Dim.Adapt & !Fast & Ready} |
             {(Alg == "DRM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "DRM") & Dim.Adapt & !Fast & !Ready} |
             {(Alg == "RAM") & !Dim.Adapt & !Fast & Ready} |
             {(Alg == "RAM") & !Dim.Adapt & !Fast & !Ready} |
             {(Alg == "RWM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "RWM") & Dim.Adapt & !Fast & !Ready}) {
               ### AMM
               if(Rec.Status > Rec.Iterations) Rec.Status <- Rec.Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AMM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, ", w=0.05))\n\n", sep="")
               }
          if({(Alg == "AM") & !Dim.Adapt & Fast & Ready} |
             {(Alg == "AM") & !Dim.Adapt & Fast & !Ready} |
             {(Alg == "AMM") & !Dim.Adapt & Fast & Ready} |
             {(Alg == "AMM") & !Dim.Adapt & Fast & !Ready} |
             {(Alg == "DRAM") & !Dim.Adapt & Fast & Ready} |
             {(Alg == "DRAM") & !Dim.Adapt & Fast & !Ready} |
             {(Alg == "MWG") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "MWG") & Dim.Adapt & !Fast & !Ready} |
             {(Alg == "RAM") & !Dim.Adapt & Fast & Ready} |
             {(Alg == "RAM") & !Dim.Adapt & Fast & !Ready}) {
               ### AMWG
               if({Alg != "AMWG"} & {Alg != "MWG"}) {
                    Rec.Status <- trunc(Rec.Status / object$Parameters)}
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"AMWG\", ",
                    "Specs=list(Periodicity=", Rec.Periodicity,
                    "))\n\n", sep="")
               }
          if((Alg == "CHARM") |
             {(Alg == "AMWG") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "AMWG") & Dim.Adapt & !Fast & !Ready} |
             {(Alg == "AMWG") & !Dim.Adapt & Fast & !Ready} | 
             {(Alg == "AMWG") & !Dim.Adapt & !Fast & !Ready} | 
             {(Alg == "HARM") & Acc.Rate.Level == 1}) {
               ### CHARM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"CHARM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if(Alg == "Slice") {
               ### Slice
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"Slice\", ",
                    "Specs=list(m=m, w=w))\n\n", sep="")
               }
          if({(Alg == "DRAM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "DRAM") & Dim.Adapt & !Fast & !Ready}) {
               ### DRAM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"DRAM\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "))\n\n", sep="")
               }
          if({(Alg == "DRM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "DRM") & Dim.Adapt & !Fast & !Ready}) {
               ### DRM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"DRM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if({(Alg == "AHMC") & Dim.Adapt & Ready} |
             {(Alg == "HMC") & Dim.Adapt & Ready}) {
               ### HMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HMC\", ",
                    "Specs=list(epsilon=", oname, "$CovarDHis[1,], ",
                    sep="")
               cat("L=", L, "))\n\n", sep="")
               }
          if(Alg == "HARM" & (Acc.Rate.Level > 1)) {
               ### HARM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HARM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if(Alg == "HMCDA" & Dim.Adapt) {
               ### HMCDA
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"HMCDA\", ",
                    "Specs=list(A=", round(Rec.Iterations/2),
                    ", delta=", delta, ",\n", sep="")
               cat("     epsilon=",
                    round(object$CovarDHis[nrow(object$CovarDHis),1],3),
                    ", Lmax=", Lmax, ", lambda=", lambda, "))\n\n", sep="")
               }
          if(Alg == "IM") {
               ### IM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"IM\", ",
                    "Specs=list(mu=apply(", oname,
                    "$Posterior1, 2, mean)))\n\n", sep="")
               }
          if(Alg == "INCA") {
               ### INCA
               library(parallel)
               detectedcores <- detectCores()
               cat(oname, " <- LaplacesDemon.hpc(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"INCA\", ",
                    "Specs=list(Adaptive=", Rec.Adaptive, ", Periodicity=",
                    Rec.Periodicity, "),\n", sep="")
               cat("     Chains=",detectedcores,", CPUs=",detectedcores,", Packages=NULL, Dyn.libs=NULL)\n\n",
                    sep="")
               }
          if({(Alg == "AMWG") & Dim.Adapt & Fast & Ready} |
             {(Alg == "AMWG") & Dim.Adapt & !Fast & Ready} |
             {(Alg == "MWG") & Dim.Adapt & Fast & Ready} |
             {(Alg == "MWG") & Dim.Adapt & !Fast & Ready}) {
               ### MWG
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"MWG\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if({Alg == "NUTS"} |
             {(Alg == "HMCDA") & !Dim.Adapt}) {
               ### NUTS
               if(Alg == "HMCDA") delta <- 0.6
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"NUTS\", ",
                    "Specs=list(A=", round(Rec.Iterations/2),
                    ", delta=", delta, ",\n", sep="")
               cat("     epsilon=",
                    round(object$CovarDHis[nrow(object$CovarDHis),1],3),
                    "))\n\n", sep="")
               }
          if({(Alg == "AM") & !Dim.Adapt & !Fast & !Ready} | 
             {(Alg == "AMM") & !Dim.Adapt & !Fast & Ready} | 
             {(Alg == "AMM") & !Dim.Adapt & !Fast & !Ready} | 
             {(Alg == "DRAM") & !Dim.Adapt & !Fast & !Ready} |
             {(Alg == "RAM") & Dim.Adapt & Fast & !Ready} |
             {(Alg == "RAM") & Dim.Adapt & !Fast & !Ready}) {
               ### RAM
               if(Rec.Status > Rec.Iterations) Rec.Status <- Rec.Iterations
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RAM\", ",
                    "Specs=list(alpha.star=0.234, Dist=\"N\", gamma=0.66,\n", sep="")
               cat("     Periodicity=10))\n\n", sep="")
               }
          if(Alg == "RJ" & (Acc.Rate.Level > 1)) {
               ### RJ
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RJ\", ",
                    "Specs=list(bin.n=bin.n, bin.p=bin.p,\n", sep="")
               cat("     parm.p=parm.p, selectable=selectable,\n", sep="")
               cat("     selected=1*(Initial.Values != 0)))\n\n", sep="")
               }
          if({(Alg == "AM") & Dim.Adapt & Fast & Ready} |
             {(Alg == "AM") & Dim.Adapt & !Fast & Ready} |
             {(Alg == "AMM") & Dim.Adapt & Fast & Ready} |
             {(Alg == "AMM") & Dim.Adapt & !Fast & Ready} |
             {(Alg == "DRAM") & Dim.Adapt & Fast & Ready} |
             {(Alg == "DRAM") & Dim.Adapt & !Fast & Ready} |
             {(Alg == "RAM") & Dim.Adapt & Fast & Ready} |
             {(Alg == "RAM") & Dim.Adapt & !Fast & Ready} |
             {(Alg == "RWM") & Dim.Adapt & Fast & Ready} |
             {(Alg == "RWM") & Dim.Adapt & !Fast & Ready}) {
               ### RWM
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"RWM\", ",
                    "Specs=NULL)\n\n", sep="")
               }
          if({(Alg == "SAMWG") & !Ready} |
             {(Alg == "SMWG") & !Ready}) {
               ### SAMWG
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"SAMWG\", ",
                    "Specs=list(Dyn=Dyn, Periodicity=", Rec.Periodicity,
                    "))\n\n", sep="")
               }
          if({(Alg == "SAMWG") & Ready} |
             {(Alg == "SMWG") & Ready}) {
               ### SMWG
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=", oname, "$Covar, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"SMWG\", ",
                    "Specs=list(Dyn=Dyn))\n\n", sep="")
               }
          if(Alg == "THMC") {
               ### THMC
               if(L > 1) {
                    L <- round(L*(Rec.Iterations/object$Iterations))
                    Rec.Iterations <- object$Iterations
                    Rec.Status <- object$Status
                    Rec.Thinning <- object$Thinning}
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"THMC\", ",
                    "Specs=list(epsilon=", oname, "$CovarDHis[1,],\n",
                    sep="")
               cat("     L=", L, ", Temperature=", Temperature,
                    "))\n\n", sep="")
               }
          if(Alg == "t-walk") {
               ### twalk
               cat(oname, " <- LaplacesDemon(Model, Data=", dname,
                    ", Initial.Values,\n", sep="")
               cat("     Covar=NULL, Iterations=",
                    Rec.Iterations, ", Status=", Rec.Status, ", ",
                    "Thinning=", Rec.Thinning, ",\n", sep="")
               cat("     Algorithm=\"twalk\", ",
                    "Specs=list(SIV=NULL, n1=4, at=6, aw=1.5))\n\n", sep="")
               }
          if((Alg == "USAMWG") | (Alg == "USMWG")) {
               ### USAMWG or USMWG
               cat("A Demonic Suggestion will not be made.\n\n")
               }
          }
     cat("Laplace's Demon is finished consorting.\n")
     }

#End

