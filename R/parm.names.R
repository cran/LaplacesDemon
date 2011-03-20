###########################################################################
# parm.names                                                              #
#                                                                         #
# The purpose of the parm.names function is to create the vector of names #
# for the parameters from a list of parameters, which may be any          #
# combination of scalars, vectors, matrices, and upper-triangular         #
# matrices. Arrays have not yet been programmed.                          #
###########################################################################

parm.names <- function(x, uppertri=NULL) {
     ### Initial Checks
     if(is.null(x)) stop("x is required in parm.names().")
     parm.length <- length(x)
     if(is.null(uppertri)) uppertri <- rep(0, parm.length)
     if(length(uppertri) != parm.length)
          stop("Length of uppertri and list attributes differs in parm.names().")
     ### Length of parm.names
     totlen <- 0
     for (i in 1:parm.length) {
          #if(is.array(x[[i]]) && (dim(x[[i]]) > 2))
          #     stop("Arrays have yet to be included in parm.names().")
          ### Scalar, Vector, or Non-Upper-Triangular Matrix
          if(uppertri[i] == 0) {
               xlen <- length(as.vector(x[[i]]))
               totlen <- totlen + xlen
               }
          ### Upper Triangular Matrix
          if(uppertri[i] == 1) {
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
         if(is.vector(x[[i]]) & (xlen > 1)) {
              for (j in 1:xlen) {
                   parm.names[cnt] <- paste(xname, "[", j, "]", sep="")
                   cnt <- cnt + 1
                   }
              }
         ### Matrix
         if(is.matrix(x[[i]]) & (uppertri[i] == 0)) {
              for (k in 1:NCOL(x[[i]])) {for (j in 1:NROW(x[[i]])) {
                   parm.names[cnt] <- paste(xname, "[", j, ",", k, "]",
                        sep="")
                   cnt <- cnt + 1
                   }}
              }
         ### Matrix, Upper Triangular
         if(is.matrix(x[[i]]) & (uppertri[i] == 1)) {
              for (k in 1:NCOL(x[[i]])) {for (j in 1:NROW(x[[i]])) {
                   if(upper.tri(x[[i]], diag=TRUE)[j,k] == TRUE) {
                        parm.names[cnt] <- paste(xname, "[", j, ",", k,
                             "]", sep="")
                        cnt <- cnt + 1
                        }
                   }}
              }
         }
     return(parm.names)
     }

### Testing
N <- 100
J <- 5
y <- rnorm(N,0,1)
X <- matrix(runif(N*J,-2,2),N,J)
S <- diag(J)
T <- diag(2)
MyData <- list(J=J, N=N, S=S, T=T, X=X, y=y,
     parm.names=parm.names(list(alpha=0, beta=rep(0,J), S=diag(J), T=diag(2)),
     uppertri=c(0,0,0,1)))
MyData
