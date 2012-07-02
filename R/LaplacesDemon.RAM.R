###########################################################################
# LaplacesDemon.RAM                                                       #
#                                                                         #
# The purpose of the LaplacesDemon.RAM function is to estimate the RAM    #
# required to update a given model and data in LaplacesDemon.             #
###########################################################################

LaplacesDemon.RAM <- function(Model, Data, Iterations, Thinning)
     {
     if(missing(Model))
          stop("The Model argument is required.")
     if(missing(Data))
          stop("The Data argument is required.")
     Const <- 1048600
     LIV <- length(Data$parm.names)
     LM <- length(Data$mon.names)
     Covar <- as.vector(object.size(matrix(runif(LIV*LIV), LIV, LIV))) /
          Const
     Data <- as.vector(object.size(Data)) / Const
     Deviance <- as.vector(object.size(runif(round(Iterations /
          Thinning)))) / Const
     Initial.Values <- as.vector(runif(LIV)) / Const
     Model <- as.vector(object.size(Model)) / Const
     Monitor <- as.vector(object.size(matrix(runif(Iterations*LM),
          round(Iterations / Thinning), LM))) / Const
     post <- as.vector(object.size(matrix(runif(Iterations*LIV),
          Iterations, LIV))) / Const
     Posterior1 <- as.vector(object.size(matrix(runif(round(Iterations /
          Thinning)), round(Iterations / Thinning), LIV))) / Const
     mem.list <- list(Covar=Covar,
          Data=Data,
          Deviance=Deviance,
          Initial.Values=Initial.Values,
          Model=Model,
          Monitor=Monitor,
          post=post,
          Posterior1=Posterior1,
          Total=sum(Covar,Data,Deviance,Initial.Values,Model,Monitor,
               post,Posterior1))
     return(mem.list)
     }

#End
