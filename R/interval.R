###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, or matrix to the interval [a,b].                      #
###########################################################################

interval <- function(x, a=-Inf, b=Inf)
     {
     if(a > b) stop("a > b in interval().")
     if(is.null(x)) stop("x is null in interval().")
     if(is.array(x) & (length(dim(x)) > 2))
          stop("arrays are unsupported by interval().")
     ### Scalar
     if(is.vector(x) && (length(x) == 1)) x <- max(a,min(b,x))
     ### Vector
     if(is.vector(x) && (length(x) > 1)) {
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)}
     ### Matrix
     if(is.matrix(x) && (length(dim(x)) == 2)) {
          r <- NROW(x)
          c <- NCOL(x)
          x <- as.vector(x)
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)
          x <- matrix(x,r,c)}
     return(x)
     }

#End
