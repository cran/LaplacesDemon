###########################################################################
# Dirichlet Distribution                                                  #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

ddirichlet <- function(x, alpha)
     {
     dirichlet1 <- function(x, alpha) {
          logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
          s <- sum((alpha-1)*log(x))
          exp(sum(s)-logD)
     }
    # make sure x is a matrix
    if(!is.matrix(x))
    if(is.data.frame(x)) x <- as.matrix(x)
    else x <- t(x)
    if(!is.matrix(alpha))
         {alpha <- matrix( alpha, ncol=length(alpha), nrow=nrow(x),
              byrow=TRUE)}
    if(any(dim(x) != dim(alpha)))
         {stop("Dimensions of x and alpha differ in ddirichlet().\n")}
    pd <- vector(length=nrow(x))
    for(i in 1:nrow(x)) {pd[i] <- dirichlet1(x[i,],alpha[i,])}
    # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
    pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] <- 0
    pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] <- 0
    return(pd)
    }

rdirichlet <- function(n, alpha)
     {
     l <- length(alpha)
     x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
     sm <- x%*%rep(1,l)
     return(x/as.vector(sm))
     }

###########################################################################
# Inverse Gamma Distribution                                              #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dinvgamma <- function(x, shape, scale = 1)
     {
     # Initial Check
     if(shape <= 0 | scale <=0) {
          stop("Shape or scale parameter negative in dinvgamma().\n")}
     alpha <- shape
     beta <- scale
     # done on log scale to allow for large alphas and betas
     log.density <- alpha * log(beta) - lgamma(alpha) -
          (alpha + 1) * log(x) - (beta/x)
     return(exp(log.density))
     }

rinvgamma <- function(n, shape, scale = 1)
     {return(1 / rgamma(n=n, shape=shape, rate=scale))}

###########################################################################
# Inverse Wishart Distribution                                            #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dinvwishart <- function(W, v, S)
     {
     if (!is.matrix(S)) S <- matrix(S)
     if (nrow(S) != ncol(S)){stop("Matrix W is not symmetric in dinvwishart().\n")}
     if (!is.matrix(W)) S <- matrix(W)
     if (nrow(W) != ncol(W)){stop("Matrix W is not symmetric in dinvwishart().\n")}
     if(nrow(S) != ncol(W)){stop("Dimensions of W and S differ in dinvwishart().\n")}
     if (v < nrow(S)){stop("v is less than the dimension of S in dinvwishart().\n")}
     k <- nrow(S)
     # denominator
     gammapart <- 1
     for(i in 1:k) {gammapart <- gammapart * gamma((v + 1 - i)/2)}
     denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
     # numerator
     detS <- det(S)
     detW <- det(W)
     hold <- S %*% solve(W)
     tracehold <- sum(hold[row(hold) == col(hold)])
     num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)
     return(num / denom)
     }

rinvwishart <- function(v, S) {return(solve(rwishart(v,solve(S))))}

###########################################################################
# Laplace Distribution                                                    #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dlaplace <- function(x, location = 0, scale = 1, log = FALSE)
     {
     logdensity = (-abs(x - location) / scale) - log(2 * scale)
     if (log == FALSE) logdensity
     else exp(logdensity)
     }
plaplace <- function(q, location = 0, scale = 1)
     {
     z = (q - location) / scale
     L = max(length(q), length(location), length(scale))
     q = rep(q, len = L)
     location = rep(location, len = L)
     scale = rep(scale, len = L)
     ifelse(q < location, 0.5 * exp(z), 1 - 0.5 * exp(-z))
     }
qlaplace <- function(p, location = 0, scale = 1)
     {
     L = max(length(p), length(location), length(scale))
     p = rep(p, len = L)
     location = rep(location, len = L)
     scale = rep(scale, len = L)
     location - sign(p - 0.5) * scale * log(2 * ifelse(p < 0.5,
          p, 1 - p))
     }
rlaplace <- function(n, location = 0, scale = 1)
     {
     location = rep(location, len = n)
     scale = rep(scale, len = n)
     r = runif(n)
     location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5,
          r, 1 - r))
     }

###########################################################################
# Multivariate Normal Distribution                                        #
#                                                                         #
# These functions are similar to those in the mnormt package.             #
###########################################################################

dmvn <- function(x, mean=rep(0,d), varcov, log=FALSE)
     {
     d  <- if(is.matrix(varcov)) ncol(varcov) else 1
     x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
     if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
     X  <- t(x - mean)
     varcov <- (varcov + t(varcov)) / 2
     u <- chol(varcov, pivot=FALSE)
     inv <- chol2inv(u)
     inv <- (inv + t(inv)) / 2
     Q  <- apply((inv %*% X)* X, 2, sum)
     logDet <- 2 * sum(log(diag(u)))
     logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
     if(log) logPDF else exp(logPDF)
     }

rmvn <- function(n=1, mean=rep(0,d), varcov)
     {
     d <- if(is.matrix(varcov)) ncol(varcov) else 1
     z <- matrix(rnorm(n*d),n,d) %*% chol(varcov)
     y <- t(mean+t(z))
     return(y)
     }

###########################################################################
# Multivariate t Distribution                                             #
#                                                                         #
# These functions are similar to those in the mnormt package.             #
###########################################################################

dmvt <- function (x, mean=rep(0,d), S, df = Inf, log = FALSE)
     {
     if (df == Inf)  return(dmvn(x, mean, S, log = log))
     d  <- if(is.matrix(S)) ncol(S) else 1
     x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
     if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
     X  <- t(x - mean)
     S <- (S + t(S)) / 2
     u <- chol(S, pivot=FALSE)
     S.inv <- chol2inv(u)
     S.inv <- (S.inv + t(S.inv)) / 2
     Q <- apply((S.inv %*% X) * X, 2, sum)
     logDet <- 2 * sum(log(diag(u)))
     logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + logDet)
          - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
     if(log) logPDF else exp(logPDF)
     }

rmvt <- function(n=1, mean=rep(0,d), S, df=Inf)
     {
     d <- if(is.matrix(S)) ncol(S) else 1
     if(df==Inf) x <- 1 else x <- rchisq(n,df)/df
     z <- rmvn(n, rep(0,d), S)
     y <- t(mean + t(z/sqrt(x)))
     return(y)
     }

###########################################################################
# Wishart Distribution                                                    #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dwishart <- function(W, v, S)
     {
     if (!is.matrix(S)) S <- matrix(S)
     if (nrow(S) != ncol(S)){stop(message="Matrix W is not symmetric in dwishart()\n\n")}
     if (!is.matrix(W)) S <- matrix(W)
     if (nrow(W) != ncol(W)){stop(message="Matrix W is not symmetric in dwishart()\n\n")}
     if(nrow(S) != ncol(W)){stop(message="Dimensions of W and S differ in dwishart()\n\n")}
     if (v < nrow(S)){stop(message="v is less than the dimension of S in dwishart()\n\n")}
     k <- nrow(S)
     # denominator
     gammapart <- 1
     for(i in 1:k) {gammapart <- gammapart * gamma((v + 1 - i)/2)}
     denom <- gammapart *  2^(v * k / 2) * pi^(k*(k-1)/4)
     # numerator
     detS <- det(S)
     detW <- det(W)
     hold <- solve(S) %*% W
     tracehold <- sum(hold[row(hold) == col(hold)])
     num <- detS^(-v/2) * detW^((v - k - 1)/2) * exp(-1/2 * tracehold)
     return(num / denom)
     }

rwishart <- function(v, S)
     {
     if (!is.matrix(S)) S <- matrix(S)
     if (nrow(S) != ncol(S)) {stop(message="Matrix S is not symmetric in rwishart().\n")}
     if (v < nrow(S)) {
          stop(message="v is less than the dimension of S in rwishart().\n")}
     p <- nrow(S)
     CC <- chol(S)
     Z <- matrix(0, p, p)
     diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
     if(p > 1)
          {
          pseq <- 1:(p-1)
          Z[rep(p*pseq, pseq) +
               unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
          }
     return(crossprod(Z %*% CC))
     }

#End
