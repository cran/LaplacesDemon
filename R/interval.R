###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, matrix, or array to the interval [a,b].               #
###########################################################################

interval <- function(x, a=-Inf, b=Inf)
     {
     if(missing(x)) stop("The x argument is required.")
     if(a > b) stop("a > b.")
     ### Scalar
     if(is.vector(x) & {length(x) == 1}) x <- max(a, min(b, x))
     ### Vector
     else if(is.vector(x) & {length(x) > 1}) {
          x.num <- which(x < a)
          x[x.num] <- a
          x.num <- which(x > b)
          x[x.num] <- b}
     ### Matrix or Array
     else if(is.array(x)) {
          d <- dim(x)
          x <- as.vector(x)
          x.num <- which(x < a)
          x[x.num] <- a
          x.num <- which(x > b)
          x[x.num] <- b
          x <- array(x, dim=d)}
     return(x)
     }

#End
