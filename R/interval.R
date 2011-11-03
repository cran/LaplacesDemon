###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, or matrix to the interval [a,b].                      #
###########################################################################

interval <- function(x, a=-Inf, b=Inf)
     {
     if(a > b) stop("a > b.")
     if(missing(x)) stop("The x argument is required.")
     if(is.array(x) & {length(dim(x)) > 2})
          stop("Arrays are unsupported.")
     ### Scalar
     if(is.vector(x) && {length(x) == 1}) x <- max(a,min(b,x))
     ### Vector
     else if(is.vector(x) && {length(x) > 1}) {
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)}
     ### Matrix
     else if(is.matrix(x) && {length(dim(x)) == 2}) {
          r <- nrow(x)
          c <- ncol(x)
          x <- as.vector(x)
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)
          x <- matrix(x,r,c)}
     ### Array (To Do)
     return(x)
     }

#End
