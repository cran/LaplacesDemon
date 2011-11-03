###########################################################################
# as.parm.names                                                           #
#                                                                         #
# The purpose of the as.parm.names function is to create the vector of    #
# names for the parameters from a list of parameters, which may be any    #
# combination of scalars, vectors, matrices, and upper-triangular         #
# matrices.                                                               #
###########################################################################

as.parm.names <- function(x, uppertri=NULL) {
     ### Initial Checks
     if(missing(x)) stop("x is required.")
     parm.length <- length(x)
     if(is.null(uppertri)) uppertri <- rep(0, parm.length)
     if(!identical(length(uppertri), parm.length))
          stop("Length of uppertri and list attributes differs.")
     ### Length of parm.names
     totlen <- 0
     for (i in 1:parm.length) {
          ### Scalar, Vector, or Non-Upper-Triangular Matrix or Aray
          if(uppertri[i] == 0) {
               xlen <- length(as.vector(x[[i]]))
               totlen <- totlen + xlen
               }
          ### Upper Triangular Matrix
          else if(uppertri[i] == 1) {
               xlen <- length(x[[i]][upper.tri(x[[i]], diag=TRUE)])
               totlen <- totlen + xlen
               }
          }
     ### Assign parm.names
     parm.names <- rep(NA, totlen)
     cnt <- 1
     for (i in 1:parm.length) {
         xname <- names(x)[i]
         xlen <- length(as.vector(x[[i]]))
         ### Scalar
         if(xlen == 1) {parm.names[cnt] <- paste(xname); cnt <- cnt + 1}
         ### Vector
         else if(is.vector(x[[i]]) & {xlen > 1}) {
              for (j in 1:xlen) {
                   parm.names[cnt] <- paste(xname, "[", j, "]", sep="")
                   cnt <- cnt + 1
                   }
              }
         ### Matrix
         else if(is.matrix(x[[i]]) & (uppertri[i] == 0)) {
              for (k in 1:ncol(x[[i]])) {for (j in 1:nrow(x[[i]])) {
                   parm.names[cnt] <- paste(xname, "[", j, ",", k, "]",
                        sep="")
                   cnt <- cnt + 1
                   }}
              }
         ### Matrix, Upper Triangular
         else if(is.matrix(x[[i]]) & (uppertri[i] == 1)) {
              for (k in 1:ncol(x[[i]])) {for (j in 1:nrow(x[[i]])) {
                   if(upper.tri(x[[i]], diag=TRUE)[j,k] == TRUE) {
                        parm.names[cnt] <- paste(xname, "[", j, ",", k,
                             "]", sep="")
                        cnt <- cnt + 1
                        }
                   }}
              }
         ### Array
         else if(is.array(x[[i]]) & (uppertri[i] == 0)) {
              count <- rep(1, length(dim(x[[i]])))
              arrayx <- array(1:prod(dim(x[[i]])), dim(x[[i]]))
              for (j in 1:prod(dim(x[[i]]))) {
                   position <- which(arrayx == j, arr.ind=TRUE)
                   parm.names[cnt] <- paste(xname, "[",
                        paste(as.vector(position), sep="", collapse=","),
                        "]", sep="")
                   cnt <- cnt + 1
                   }
              }
         ### Array, Upper Triangular
         else if(is.array(x[[i]]) & (uppertri[i] == 1)) {
              stop("upper.tri does not function with arrays.")}
         }
     return(parm.names)
     }

#End
