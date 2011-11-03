###########################################################################
# BayesTheorem                                                            #
###########################################################################

BayesTheorem <- function(PrA, PrBA)
     {
     if(missing(PrA)) stop("The PrA argument is required.")
     if(missing(PrBA)) stop("The PrBA argument is required.")
     PrAB <- (PrBA * PrA) / sum(PrBA * PrA)
     class(PrAB) <- "bayestheorem"
     return(PrAB)
     }

#End
