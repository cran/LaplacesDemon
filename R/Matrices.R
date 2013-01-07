###########################################################################
# Matrices                                                                #
#                                                                         #
# These are utility functions for matrices.                               #
###########################################################################

as.indicator.matrix <- function(x)
     {
     n <- length(x)
     x <- as.factor(x)
     X <- matrix(0, n, length(levels(x)))
     X[(1:n) + n*(unclass(x)-1)] <- 1
     dimnames(X) <- list(names(x), levels(x))
     return(X)
     }
as.inverse <- function(x)
     {
     if(!is.matrix(x)) x <- matrix(x)
     if(!is.square.matrix(x)) stop("x must be a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x must be a symmetric matrix.")
     tol <- .Machine$double.eps
     options(show.error.messages=FALSE)
     xinv <- try(solve(x))
     if(inherits(xinv, "try-error")) {
          k <- nrow(x)
          eigs <- eigen(x, symmetric=TRUE)
          if(min(eigs$values) < tol) {
               tolmat <- diag(k)
               for (i in 1:k)
                    if(eigs$values[i] < tol) tolmat[i,i] <- 1/tol
                    else tolmat[i,i] <- 1/eigs$values[i]
               }
          else tolmat <- diag(1/eigs$values, nrow=length(eigs$values))
          xinv <- eigs$vectors %*% tolmat %*% t(eigs$vectors)
          }
     options(show.error.messages=TRUE)
     xinv <- as.symmetric.matrix(xinv)
     return(xinv)
     }
as.parm.matrix <- function(x, k, parm, Data, a=-Inf, b=Inf, restrict=FALSE,
     chol=FALSE)
     {
     X <- matrix(0, k, k)
     if(restrict == TRUE) {
          X[upper.tri(X, diag=TRUE)] <- c(1,
               parm[grep(deparse(substitute(x)),
               Data$parm.names)])}
     else {
          X[upper.tri(X, diag=TRUE)] <- parm[grep(deparse(substitute(x)),
               Data$parm.names)]}
     if(chol == TRUE) {
          if(a != -Inf | b != Inf) {
               x <- as.vector(X[upper.tri(X, diag=TRUE)])
               x.num <- which(x < a)
               x[x.num] <- a
               x.num <- which(x > b)
               x[x.num] <- b
               X[upper.tri(X, diag=TRUE)] <- x
               diag(X) <- abs(diag(X))
               }
          X[lower.tri(X)] <- 0
          return(X)
          }
     X[lower.tri(X)] <- t(X)[lower.tri(X)]
     if(a != -Inf | b != Inf) {
          x <- as.vector(X[upper.tri(X, diag=TRUE)])
          x.num <- which(x < a)
          x[x.num] <- a
          x.num <- which(x > b)
          x[x.num] <- b
          X[upper.tri(X, diag=TRUE)] <- x
          X[lower.tri(X)] <- t(X)[lower.tri(X)]
          }
     if(!is.symmetric.matrix(X)) X <- as.symmetric.matrix(X)
     if(!exists("LDEnv")) LDEnv <- new.env()
     if(restrict == FALSE) {
          if(is.positive.definite(X)) {
               assign("LaplacesDemonMatrix", as.vector(X[upper.tri(X,
                    diag=TRUE)]), envir=LDEnv)}
          else {
               if(exists("LaplacesDemonMatrix", envir=LDEnv)) {
                    X[upper.tri(X,
                         diag=TRUE)] <- as.vector(get("LaplacesDemonMatrix",
                         envir=LDEnv))
                    X[lower.tri(X)] <- t(X)[lower.tri(X)]}
               else {X <- diag(k)}}
          }
     if(restrict == TRUE) {
          if(is.positive.definite(X)) {
               assign("LaplacesDemonMatrix", as.vector(X[upper.tri(X,
                    diag=TRUE)][-1]), envir=LDEnv)}
          else {
               if(exists("LaplacesDemonMatrix", envir=LDEnv)) {
                    X[upper.tri(X, diag=TRUE)] <- c(1,
                         as.vector(get("LaplacesDemonMatrix",
                         envir=LDEnv)))
                    X[lower.tri(X)] <- t(X)[lower.tri(X)]
                    if(!is.symmetric.matrix(X)) X <- as.symmetric.matrix(X)
                    }
               else {X <- diag(k)}}
          }
     return(X)
     }
as.positive.definite <- function(x)
     {
     eig.tol <- 1e-06
     conv.tol <- 1e-07
     posd.tol <- 1e-08
     iter <- 0; maxit <- 100
     n <- ncol(x)
     D_S <- x
     D_S[] <- 0
     X <- x
     converged <- FALSE
     conv <- Inf
     while (iter < maxit && !converged) {
          Y <- X
          R <- Y - D_S
          e <- eigen(R, symmetric=TRUE)
          Q <- e$vectors
          d <- e$values
          p <- d > eig.tol * d[1]
          if(!any(p))
               stop("Matrix seems negative semi-definite.")
          Q <- Q[, p, drop=FALSE]
          X <- tcrossprod(Q * rep(d[p], each=nrow(Q)), Q)
          D_S <- X - R
          conv <- norm(Y - X, "I") / norm(Y, "I")
          iter <- iter + 1
          converged <- (conv <= conv.tol)
          }
     if(!converged) {
          warning("as.positive.definite did not converge in ", iter,
               " iterations.")}
     e <- eigen(X, symmetric=TRUE)
     d <- e$values
     Eps <- posd.tol * abs(d[1])
     if(d[n] < Eps) {
          d[d < Eps] <- Eps
          Q <- e$vectors
          o.diag <- diag(X)
          X <- Q %*% (d * t(Q))
          D <- sqrt(pmax(Eps, o.diag)/diag(X))
          X[] <- D * X * rep(D, each=n)}
     X <- as.symmetric.matrix(X)
     return(X)
     }
as.positive.semidefinite <- function(x)
     {
     if(!is.matrix(x)) x <- matrix(x)
     if(!is.square.matrix(x)) stop("x must be a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x must be a symmetric matrix.")
     iter <- 0; maxit <- 100
     converged <- FALSE
     while (iter < maxit && !converged) {
          iter <- iter + 1
          out <- eigen(x=x, symmetric=TRUE)
          mGamma <- t(out$vectors)
          vLambda <- out$values
          vLambda[vLambda < 0] <- 0
          x <- t(mGamma) %*% diag(vLambda) %*% mGamma
          x <- as.symmetric.matrix(x)
          if(is.positive.semidefinite(x)) converged <- TRUE
          }
     if(converged == FALSE) {
          warning("as.positive.semidefinite did not converge in ", iter,
               " iterations.")}
     return(x)
     }
as.symmetric.matrix <- function(x, k=NULL)
     {
     if(is.vector(x)) {
          if(any(!is.finite(x))) stop("x must have finite values.")
          if(is.null(k)) k <- (-1 + sqrt(1 + 8 * length(x))) / 2
          symm <- matrix(0, k, k)
          symm[lower.tri(symm, diag=TRUE)] <- x
          symm2 <- symm
          symm2[upper.tri(symm2, diag=TRUE)] <- 0
          symm <- symm + t(symm2)
          }
     else if(is.matrix(x)) {
          if(!is.square.matrix(x)) stop("x must be a square matrix.")
          if(any(!is.finite(diag(x))))
               stop("The diagonal of x must have finite values.")
          symm <- x
          x.lower.fin <- FALSE; x.upper.fin <- FALSE
          if(all(is.finite(x[lower.tri(x, diag=TRUE)])))
               x.lower.fin <- TRUE
          if(all(is.finite(x[upper.tri(x, diag=TRUE)])))
               x.upper.fin <- TRUE
          if(x.lower.fin) symm[upper.tri(x)] <- t(x)[upper.tri(x)]
          else if(x.upper.fin) symm[lower.tri(x)] <- t(x)[lower.tri(x)]
          else {
               new.up <- ifelse(!is.finite(x[upper.tri(x)]),
                    t(x)[lower.tri(x)], x[upper.tri(x)])
               new.low <- ifelse(!is.finite(x[lower.tri(x)]),
                    t(x)[upper.tri(x)], x[lower.tri(x)])
               if(any(!is.finite(c(new.up, new.low))))
                    stop("Off-diagonals in x must have finite values.")
               else {
                    symm[upper.tri(symm)] <- new.up
                    symm[lower.tri(symm)] <- new.low
                    }
               }
     }
     else stop("x must be a vector or matrix.")
     return(symm)
     }
Cov2Cor <- function(Sigma)
     {
     if(missing(Sigma)) stop("Sigma is a required argument.")
     if(any(!is.finite(Sigma))) stop("Sigma must have finite values.")
     if(is.matrix(Sigma)) {
          if(!is.positive.definite(Sigma))
               stop("Sigma is not positive-definite.")
          x <- 1 / sqrt(diag(Sigma))
          R <- x * t(x * Sigma)}
     else if(is.vector(Sigma)) {
          k <- as.integer(sqrt(length(Sigma)))
          Sigma <- matrix(Sigma, k, k)
          x <- 1 / sqrt(diag(Sigma))
          R <- as.vector(x * t(x * Sigma))}
     return(R)
     }
Hessian <- function(Model, parm, Data)
     {
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(parm)) stop("The parm argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     Interval <- 1.0E-6
     parm.len <- length(parm)
     eps <- Interval * parm
     H <- matrix(0, parm.len, parm.len)
     for (i in 1:parm.len) {
          for (j in i:parm.len) {
               x1 <- x2 <- x3 <- x4 <- parm
               x1[i] <- x1[i] + eps[i]
               x1[j] <- x1[j] + eps[j]
               x2[i] <- x2[i] + eps[i]
               x2[j] <- x2[j] - eps[j]
               x3[i] <- x3[i] - eps[i]
               x3[j] <- x3[j] + eps[j]
               x4[i] <- x4[i] - eps[i]
               x4[j] <- x4[j] - eps[j]
               H[i, j] <- {Model(x1, Data)[[1]] -
                    Model(x2, Data)[[1]] - Model(x3, Data)[[1]] +
                    Model(x4, Data)[[1]]} / {4 * eps[i] * eps[j]}
               }
          }
     H[lower.tri(H)] <- t(H)[lower.tri(H)]
     return(H)
     }
is.positive.definite <- function(x)
     {
     if(!is.matrix(x)) stop("x is not a matrix.")
     if(!is.square.matrix(x)) stop("x is not a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x is not a symmetric matrix.")
     ### Deprecated Method 1
     #pd <- TRUE
     #ed <- eigen(x, symmetric=TRUE)
     #ev <- ed$values
     #if(!all(ev >= -1e-06 * abs(ev[1]))) pd <- FALSE
     ### Deprecated Method 2
     #eval <- eigen(x, only.values=TRUE)$values
     #if(any(is.complex(eval))) eval <- rep(0, max(dim(x)))
     #tol <- max(dim(x)) * max(abs(eval)) * .Machine$double.eps
     #if(all(eval > tol)) pd <- TRUE
     ### Currently Active, Method 3
     eigs <- eigen(x, symmetric=TRUE)$values
     if(any(is.complex(eigs))) return(FALSE)
     if(all(eigs > 0)) pd <- TRUE
     else pd <- FALSE
     return(pd)
     }
is.positive.semidefinite <- function(x)
     {
     if(!is.matrix(x)) stop("x is not a matrix.")
     if(!is.square.matrix(x)) stop("x is not a square matrix.")
     if(!is.symmetric.matrix(x)) stop("x is not a symmetric matrix.")
     eigs <- eigen(x, symmetric=TRUE)$values
     if(any(is.complex(eigs))) return(FALSE)
     if(all(eigs >= 0)) pd <- TRUE
     else pd <- FALSE
     return(pd)
     }
is.square.matrix <- function(x) {return(nrow(x) == ncol(x))}
is.symmetric.matrix <- function(x) {return(sum(x == t(x)) == (nrow(x)^2))}
lower.triangle <- function(x, diag=FALSE)
     {
     return(x[lower.tri(x, diag=diag)])
     }
read.matrix <- function(file, header=FALSE, sep="", skip=0, na.rm=FALSE) 
     {
     row.lens <- count.fields(file, sep=sep, skip=skip)
     if(any(row.lens != row.lens[1])) 
          stop("The number of columns is not constant.")
     if(header) {
          nrows <- length(row.lens) - 1
          ncols <- row.lens[2]
          col.names <- scan(file, what="", sep=sep, nlines=1, 
            quiet=TRUE, skip=skip)
          x <- scan(file, sep=sep, skip=skip + 1, quiet=TRUE)}
     else {
          nrows <- length(row.lens)
          ncols <- row.lens[1]
          x <- scan(file, sep=sep, skip=skip, quiet=TRUE)
          col.names <- NULL}
     x <- as.double(x)
     if(ncols > 1) {
          dim(x) <- c(ncols, nrows)
          x <- t(x)
          colnames(x) <- col.names}
     else if(ncols == 1) x <- as.vector(x)
     else stop("The number of columns is wrong.")
     num.mis <- sum(is.na(x))
     if(num.mis > 0) {
          cat("\n", num.mis, "missing value(s) found.")
          cat("\n", sum(complete.cases(x)),
               "row(s) found with missing values.")}
     if(na.rm == TRUE) x <- x[complete.cases(x),]
     return(x)
     }
tr <- function(x) {return(sum(diag(x)))}
upper.triangle <- function(x, diag=FALSE)
     {
     return(x[upper.tri(x, diag=diag)])
     }

#End
