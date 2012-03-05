###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, matrix, or array to the interval [a,b].               #
###########################################################################

interval <- function(x, a=-Inf, b=Inf)
     {
     if(a > b) stop("a > b.")
     if(missing(x)) stop("The x argument is required.")
     ### Scalar
     if(is.vector(x) && {length(x) == 1}) x <- max(a,min(b,x))
     ### Vector
     else if(is.vector(x) && {length(x) > 1}) {
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)}
     ### Matrix or Array
     else if(is.array(x)) {
          d <- dim(x)
          x <- as.vector(x)
          x <- ifelse(x < a, a, x)
          x <- ifelse(x > b, b, x)
          x <- array(x, dim=d)}
     return(x)
     }

#End
