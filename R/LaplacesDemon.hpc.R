###########################################################################
# LaplacesDemon.hpc                                                       #
#                                                                         #
# The purpose of the LaplacesDemon.hpc function is to extend              #
# LaplacesDemon to parallel processing on multiple cores.                 #
###########################################################################

#I don't know how to 'shut up' the makeCluster function...
#no luck with capture.output, sink, or invisible

#I'd like to send different sets of initial values to different clusters.

LaplacesDemon.hpc <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=100000, Status=1000, Thinning=100, Algorithm="RWM",
     Specs=NULL, Chains=2, CPUs=2, Packages=NULL, Dyn.libs=NULL)
     {
     library(parallel, quietly=TRUE)
     detectedCores <- detectCores()
     cat("\n\nCPUs Detected:", detectedCores,"\n")
     if(CPUs > detectedCores) {
          cat("\nOnly", detectedCores, "will be used.\n")
          CPUs <- detectedCores}
     cat("\nLaplace's Demon is preparing environments for CPUs...")
     cat("\n##################################################\n")
     cl <- makeCluster(detectedCores)
     cat("\n##################################################\n")
     on.exit({stopCluster(cl); cat("\n\nLaplace's Demon has finished.\n")})
     varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
          ls(envir=parent.env(environment()))))
     clusterExport(cl, varlist=varlist, envir=environment())
     clusterSetRNGStream(cl)
     wd <- getwd()
     clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
          envir=environment())
     demon.wrapper <- function(x, ...)
          {
          if(!is.null(Packages)) {
               sapply(Packages,
                    function(x) library(x, character.only=TRUE,
                         quietly=TRUE))}
          if(!is.null(Dyn.libs)) {
               sapply(Dyn.libs,
                    function(x) dyn.load(paste(wd, x, sep = "/")))
               on.exit(sapply(Dyn.libs,
                    function(x) dyn.unload(paste(wd, x, sep = "/"))))}
          LaplacesDemon(...)
          }
     cat("\nStatus messages are not displayed for parallel processing.")
     cat("\nLaplace's Demon is beginning parallelization...\n")
     LaplacesDemon.out <- clusterApply(cl, 1:Chains, demon.wrapper, Model,
          Data, Initial.Values, Covar, Iterations, Status, Thinning,
          Algorithm, Specs)
     class(LaplacesDemon.out) <- "demonoid.hpc"
     return(LaplacesDemon.out)
     }

#End
