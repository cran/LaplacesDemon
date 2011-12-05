###########################################################################
# GIV                                                                     #
#                                                                         #
# GIV stands for ``generate initial values'', and the purpose of the GIV  #
# function is to generate initial values for the LaplaceApproximation and #
# LaplacesDemon functions.                                                #
###########################################################################

GIV <- function(Model, Data, n=1000)
     {
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data$parm.names)) stop("parm.names missing in Data.")
     LIV <- length(Data$parm.names)
     high <- 100; low <- -100
     a <- Model(rep(low, LIV), Data)[[5]]
     b <- Model(rep(high, LIV), Data)[[5]]
     ab.range <- b - a
     ab.mu <- a + ab.range / 2
     ab.mu <- ifelse({a == 0} & {b == high}, 10, ab.mu)
     Scale <- 1 / ab.range
     Scale <- ifelse(ab.range == 0, 0, Scale)
     iv <- rep(NA, LIV)
     for (i in 1:n) {
          IV <- rnorm(LIV, ab.mu, ab.range * Scale)
          M <- Model(IV, Data)
          if(is.finite(M[[1]]) & is.finite(M[[2]])) {
               iv <- IV; break}
          Scale <- Scale + ab.range / n / 2}
     if(any(is.na(iv))) 
          cat("\nWARNING: Acceptable initial values were not generated.")
     return(iv)
     }
#End
