###########################################################################
# Dirichlet Distribution                                                  #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

ddirichlet <- function(x, alpha, log=FALSE)
     {
     dirichlet1 <- function(x, alpha) {
          logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
          s <- sum((alpha-1)*log(x))
          exp(sum(s)-logD)
          }
    if(!is.matrix(x)) x <- matrix(x)
    if(is.data.frame(x)) x <- matrix(x)
    else x <- t(x)
    if(!is.matrix(alpha))
         {alpha <- matrix(alpha, ncol=length(alpha), nrow=nrow(x),
              byrow=TRUE)}
    if(any(dim(x) != dim(alpha)))
         {stop("Dimensions of x and alpha differ in ddirichlet().\n")}
    dens <- vector(length=nrow(x))
    for (i in 1:nrow(x)) {dens[i] <- dirichlet1(x[i,], alpha[i,])}
    # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
    dens[apply(x, 1, function(z) any(z <0 | z > 1))] <- 0
    dens[apply(x, 1, function(z) all.equal(sum(z), 1) !=TRUE)] <- 0
    if(log == TRUE) dens <- log(dens)
    return(dens)
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

dinvgamma <- function(x, shape, scale=1, log=FALSE)
     {
     # Initial Check
     if(shape <= 0 | scale <=0) {
          stop("Shape or scale parameter negative in dinvgamma().\n")}
     alpha <- shape
     beta <- scale
     # done on log scale to allow for large alphas and betas
     log.density <- alpha * log(beta) - lgamma(alpha) -
          (alpha + 1) * log(x) - (beta/x)
     dens <- exp(log.density)
     if(log == TRUE) dens <- log.density
     return(dens)
     }

rinvgamma <- function(n, shape, scale=1)
     {return(1 / rgamma(n=n, shape=shape, rate=scale))}

###########################################################################
# Inverse Wishart Distribution                                            #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dinvwishart <- function(Sigma, nu, R, log=FALSE)
     {
     if(!is.matrix(R)) R <- matrix(R)
     if(nrow(R) != ncol(R)) {stop("Matrix R is not symmetric in dinvwishart().\n")}
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(nrow(Sigma) != ncol(Sigma)) {stop("Matrix Sigma is not symmetric in dinvwishart().\n")}
     if(nrow(R) != ncol(Sigma)) {stop("Dimensions of Sigma and R differ in dinvwishart().\n")}
     if(nu < nrow(R)) {stop("nu is less than the dimension of R in dinvwishart().\n")}
     p <- nrow(R)
     # denominator
     gammapart <- 1
     for (i in 1:p) {gammapart <- gammapart * gamma((nu + 1 - i)/2)}
     denom <- gammapart *  2^(nu * p / 2) * pi^(p*(p-1)/4)
     # numerator
     detR <- det(R)
     detSigma <- det(Sigma)
     hold <- R %*% solve(Sigma)
     tracehold <- sum(hold[row(hold) == col(hold)])
     num <- detR^(nu/2) * detSigma^(-(nu + p + 1)/2) * exp(-1/2 * tracehold)
     dens <- num / denom
     if(log == TRUE) dens <- log(num / denom)
     return(dens)
     }

rinvwishart <- function(nu, R) {return(solve(rwishart(nu,solve(R))))}

###########################################################################
# Laplace Distribution                                                    #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dlaplace <- function(x, location=0, scale=1, log=FALSE)
     {
     if(scale <= 0) {stop("scale parameter negative in dlaplace().\n")}
     logdens = (-abs(x - location) / scale) - log(2 * scale)
     if(log == TRUE) logdens
     else exp(logdens)
     }
plaplace <- function(q, location=0, scale=1)
     {
     if(scale <= 0) {stop("scale parameter negative in plaplace().\n")}
     z = (q - location) / scale
     L = max(length(q), length(location), length(scale))
     q = rep(q, len=L)
     location = rep(location, len=L)
     scale = rep(scale, len=L)
     ifelse(q < location, 0.5 * exp(z), 1 - 0.5 * exp(-z))
     }
qlaplace <- function(p, location=0, scale=1)
     {
     if(scale <= 0) {stop("scale parameter negative in qlaplace().\n")}
     L = max(length(p), length(location), length(scale))
     p = rep(p, len=L)
     location = rep(location, len=L)
     scale = rep(scale, len=L)
     location - sign(p - 0.5) * scale * log(2 * ifelse(p < 0.5,
          p, 1 - p))
     }
rlaplace <- function(n, location=0, scale=1)
     {
     if(scale <= 0) {stop("scale parameter negative in rlaplace().\n")}
     location = rep(location, len=n)
     scale = rep(scale, len=n)
     r = runif(n)
     location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5,
          r, 1 - r))
     }

###########################################################################
# Multivariate Normal Distribution                                        #
#                                                                         #
# These functions are similar to those in the mnormt and mvtnorm          #
# packages.                                                               #
###########################################################################

dmvn <- function(x, mu=rep(0,d), Sigma, log=FALSE)
     {
     ### Like in the mnormt package...
     #d  <- if(is.matrix(Sigma)) ncol(Sigma) else 1
     #x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
     #if(is.vector(mu)) mu <- outer(rep(1, nrow(x)), mu)
     #X  <- t(x - mu)
     #Sigma <- (Sigma + t(Sigma)) / 2
     #u <- chol(Sigma, pivot=FALSE)
     #inv <- chol2inv(u)
     #inv <- (inv + t(inv)) / 2
     #Q  <- apply((inv %*% X)* X, 2, sum)
     #logDet <- 2 * sum(log(diag(u)))
     #logdens <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
     #if(log) logdens else exp(logdens)
     ### Like in the mvtnorm package...
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     d  <- if(is.matrix(Sigma)) ncol(Sigma) else 1
     options(warn=-1)
     distval <- mahalanobis(x, center = mu, cov = Sigma)
     logdet <- sum(log(eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values))
     logdens <- -(ncol(x)*log(2*pi) + logdet + distval)/2
     options(warn=0)
     if(log) return(logdens)
     exp(logdens)
     }

rmvn <- function(n=1, mu=rep(0,d), Sigma)
     {
     d <- if(is.matrix(Sigma)) ncol(Sigma) else 1
     z <- matrix(rnorm(n*d),n,d) %*% chol(Sigma)
     y <- t(mu + t(z))
     return(y)
     }

###########################################################################
# Multivariate t Distribution                                             #
#                                                                         #
# These functions are similar to those in the mnormt package.             #
###########################################################################

dmvt <- function (x, mu=rep(0,d), S, df=Inf, log=FALSE)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(df == Inf) return(dmvn(x, mu, S, log = log))
     d  <- if(is.matrix(S)) ncol(S) else 1
     x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
     if(is.vector(mu)) mu <- outer(rep(1, nrow(x)), mu)
     X  <- t(x - mu)
     S <- (S + t(S)) / 2
     u <- chol(S, pivot=FALSE)
     S.inv <- chol2inv(u)
     S.inv <- (S.inv + t(S.inv)) / 2
     Q <- apply((S.inv %*% X) * X, 2, sum)
     logDet <- 2 * sum(log(diag(u)))
     logdens <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + logDet)
          - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
     if(log) logdens else exp(logdens)
     }

rmvt <- function(n=1, mu=rep(0,d), S, df=Inf)
     {
     if(!is.matrix(S)) S <- matrix(S)
     d <- if(is.matrix(S)) ncol(S) else 1
     if(df==Inf) x <- 1 else x <- rchisq(n,df)/df
     z <- rmvn(n, rep(0,d), S)
     y <- t(mu + t(z/sqrt(x)))
     return(y)
     }

###########################################################################
# Wishart Distribution                                                    #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dwishart <- function(Omega, nu, S, log=FALSE)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(nrow(S) != ncol(S)) {stop("Matrix Omega is not symmetric in dwishart()\n\n")}
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(nrow(Omega) != ncol(Omega)) {stop("Matrix Omega is not symmetric in dwishart()\n\n")}
     if(nrow(S) != ncol(Omega)) {stop("Dimensions of Omega and S differ in dwishart()\n\n")}
     if(nu < nrow(S)) {stop("nu is less than the dimension of S in dwishart()\n\n")}
     p <- nrow(S)
     # denominator
     gammapart <- 1
     for (i in 1:p) {gammapart <- gammapart * gamma((nu + 1 - i)/2)}
     denom <- gammapart * 2^(nu * p / 2) * pi^(p*(p-1)/4)
     # numerator
     detS <- det(S)
     detOmega <- det(Omega)
     hold <- solve(S) %*% Omega
     tracehold <- sum(hold[row(hold) == col(hold)])
     num <- detS^(-nu/2) * detOmega^((nu - p - 1)/2) * exp(-1/2 * tracehold)
     dens <- num / denom
     if (log == TRUE) dens <- log(num / denom)
     return(dens)
     }

rwishart <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(nrow(S) != ncol(S)) {stop("Matrix S is not symmetric in rwishart().\n")}
     if(nu < nrow(S)) {
          stop("nu is less than the dimension of S in rwishart().\n")}
     p <- nrow(S)
     CC <- chol(S)
     Z <- matrix(0, p, p)
     diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
     if(p > 1)
          {
          pseq <- 1:(p-1)
          Z[rep(p*pseq, pseq) +
               unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
          }
     return(crossprod(Z %*% CC))
     }

#End
