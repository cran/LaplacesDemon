###########################################################################
# Cov2Cor                                                                 #
#                                                                         #
# The purpose of the Cov2Cor function is to convert a k x k covariance    #
# matrix to a correlation matrix. The covariance matrix may be supplied   #
# either in matrix or vector form, and the correlation matrix will be     #
# returned accordingly in matrix or vector form.                          #
###########################################################################

Cov2Cor <- function(Sigma)
     {
     if(is.matrix(Sigma)) {
          if(nrow(Sigma) != ncol(Sigma)) stop("Sigma is not symmetric in cov2cor().")
          x <- 1 / sqrt(diag(Sigma))
          R <- x * t(x * Sigma)}
     if(is.vector(Sigma)) {
          k <- as.integer(sqrt(length(Sigma)))
          Sigma <- matrix(Sigma, k, k)
          x <- 1 / sqrt(diag(Sigma))
          R <- as.vector(x * t(x * Sigma))}
     return(R)
     }

#End

