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
     if(any(is.na(Sigma))) stop("Sigma has missing values in Cov2Cor().")
     if(any(is.nan(Sigma)))
          stop("Sigma has non-numeric values (NaN's) in Cov2Cor().")
     if(any(is.infinite(Sigma)))
          stop("Sigma has infinite values in Cov2Cor().")
     if(is.matrix(Sigma)) {
          if(nrow(Sigma) != ncol(Sigma))
               stop("Sigma is not symmetric in Cov2Cor().")
          x <- 1 / sqrt(diag(Sigma))
          R <- x * t(x * Sigma)}
     else if(is.vector(Sigma)) {
          k <- as.integer(sqrt(length(Sigma)))
          Sigma <- matrix(Sigma, k, k)
          x <- 1 / sqrt(diag(Sigma))
          R <- as.vector(x * t(x * Sigma))}
     return(R)
     }

#End

