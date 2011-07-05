###########################################################################
# Bernoulli Distribution                                                  #
#                                                                         #
# These functions are similar to those in the Rlab package.               #
###########################################################################

dbern <- function(x, prob, log=FALSE)
     {dbinom(x, 1, prob, log)}
pbern <- function(q, prob, lower.tail=TRUE, log.p=FALSE)
     {pbinom(q, 1, prob, lower.tail, log.p)}
qbern <- function(p, prob, lower.tail=TRUE, log.p=FALSE)
     {qbinom(p, 1, prob, lower.tail, log.p)}
rbern <- function(n, prob)
     {rbinom(n, 1, prob)}

###########################################################################
# Categorical Distribution                                                #
#                                                                         #
###########################################################################

dcat <- function(x, p, log=FALSE)
     {
     if(!is.matrix(p)) p <- rbind(p)
     if(is.vector(x) & (length(x) == 1)) {
          temp <- rep(0, ncol(p))
          temp[x] <- 1
          x <- t(temp)
          }
     if(is.vector(x) & (length(x) > 1)) x <- indmat(x)
     if(!identical(dim(x),dim(p))) stop("Dimensions of x and p differ in dcat().")
     dens <- x*p
     if(log == TRUE) dens <- x*log(p)
     dens <- as.vector(rowSums(dens))
     return(dens)
     }

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
     x <- matrix(rgamma(l*n,alpha), ncol=l, byrow=TRUE)
     sm <- x %*% rep(1,l)
     return(x/as.vector(sm))
     }

###########################################################################
# Half-Cauchy Distribution                                                #
#                                                                         #
###########################################################################

dhalfcauchy <- function(x, scale=25, log=FALSE)
     {
     x <- as.vector(x)
     if(scale <= 0) {stop("scale parameter negative in dhalfcauchy().\n")}
     dens = 2*scale / (pi*(x^2 + scale^2))
     if(log == TRUE) dens <- log(dens)
     return(dens)
     }
phalfcauchy <- function(q, scale=25)
     {
     q <- as.vector(q)
     if(scale <= 0) {stop("scale parameter negative in phalfcauchy().\n")}
     z <- (2/pi)*atan(q/scale)
     return(z)
     }
qhalfcauchy <- function(p, scale=25)
     {
     p <- as.vector(p)
     if(scale <= 0) {stop("scale parameter negative in qhalfcauchy().\n")}
     x <- scale*tan((pi*p)/2)
     return(x)
     }
rhalfcauchy <- function(n, scale=25)
     {
     n <- as.vector(n)
     if(scale <= 0) {stop("scale parameter negative in rhalfcauchy().\n")}
     p <- runif(n, 0, 1)
     x <- scale*tan((pi*p)/2)
     return(x)
     }

###########################################################################
# Half-t Distribution                                                     #
#                                                                         #
###########################################################################

dhalft <- function(x, scale=25, nu=1, log=FALSE)
     {
     x <- as.vector(x)
     if(scale <= 0) {stop("scale parameter negative in dhalft().\n")}
     dens = (1 + (1/nu)*(x/scale)^2)^(-(nu+1)/2)
     if(log == TRUE) dens <- log(dens)
     return(dens)
     }
phalft <- function(q, scale=25, nu=1)
     {
     q <- as.vector(q)
     if(scale <= 0) {stop("scale parameter negative in phalft().\n")}
     z <- ptrunc(q, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(z)
     }
qhalft <- function(p, scale=25, nu=1)
     {
     p <- as.vector(p)
     if(scale <= 0) {stop("scale parameter negative in qhalft().\n")}
     x <- rtrunc(p, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(x)
     }
rhalft <- function(n, scale=25, nu=1)
     {
     n <- as.vector(n)
     if(scale <= 0) {stop("scale parameter negative in rhalft().\n")}
     x <- rtrunc(n, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(x)
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
     x <- as.vector(x)
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
     x <- as.vector(x); location <- as.vector(location)
     if(scale <= 0) {stop("scale parameter negative in dlaplace().\n")}
     logdens = (-abs(x - location) / scale) - log(2 * scale)
     if(log == FALSE) logdens <- exp(logdens)
     }
plaplace <- function(q, location=0, scale=1)
     {
     q <- as.vector(q); location <- as.vector(location)
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
     p <- as.vector(p); location <- as.vector(location)
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

dmvn <- function(x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     p <- nrow(Sigma)
     ed <- eigen(Sigma, symmetric = TRUE)
     ev <- ed$values
     #if(!all(ev >= -1e-06 * abs(ev[1])))
     #     stop("Sigma is not positive definite in dmvn().")
     ss <- if(!is.matrix(mu)) {
          x - rep(mu, each = nrow(x))
          } else {
               x - mu}
     inv.Sigma <- ed$vectors %*% (t(ed$vectors)/ev)
     quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
     options(warn=-1)
     if(is.nan(sum(log(ev)))) z <- log(1.0E-100) else z <- sum(log(ev))
     options(warn=0)
     fact <- -0.5 * (p * log(2 * pi) + z)
     dens <- as.vector(fact - quad)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }

rmvn <- function(n=1, mu=rep(0,k), Sigma)
     {
     k <- if(is.matrix(Sigma)) ncol(Sigma) else 1
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     y <- t(mu + t(z))
     return(y)
     }

###########################################################################
# Multivariate t Distribution                                             #
#                                                                         #
# These functions are similar to those in the mnormt package.             #
###########################################################################

dmvt <- function(x, mu, S, df=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     p <- nrow(S)
     ed <- eigen(S, symmetric = TRUE)
     ev <- ed$values
     #if (!all(ev >= -1e-06 * abs(ev[1])))
     #    stop("S is not positive definite in dmvt().")
     ss <- if(!is.matrix(mu)) {
          x - rep(mu, each = nrow(x))
          } else {
               x - mu}
     inv.S <- ed$vectors %*% (t(ed$vectors)/ev)
     quad <- rowSums((ss %*% inv.S) * ss)/df
     options(warn=-1)
     if(is.nan(sum(log(ev)))) z <- log(1.0E-100) else z <- sum(log(ev))
     options(warn=0)
     fact <- lgamma((df + p)/2) - lgamma(df/2) -
          0.5 * (p * (log(pi) + log(df)) + z)
     dens <- fact - 0.5 * (df + p) * log(1 + quad)
     if(log == FALSE) dens <- exp(fact) * ((1 + quad)^(-(df + p)/2))
     return(dens)
     }

rmvt <- function(n=1, mu=rep(0,k), S, df=Inf)
     {
     if(!is.matrix(S)) S <- matrix(S)
     k <- if(is.matrix(S)) ncol(S) else 1
     if(df==Inf) x <- 1 else x <- rchisq(n,df)/df
     z <- rmvn(n, rep(0,k), S)
     y <- t(mu + t(z/sqrt(x)))
     return(y)
     }

###########################################################################
# Student t distribution (3-parameter)                                    #
#                                                                         #
# These functions are similar to those in the gamlss.dist package.        #
###########################################################################

dst<-function(x, mu=0, sigma=1, nu=10, log=FALSE)
     {
     if(any(sigma <= 0)) stop("sigma must be positive in dst().")
     if(any(nu <= 0))  stop("nu must be positive in dst().")
     if(length(nu) > 1) dens <- ifelse(nu > 1000000,
          dnorm(x, mu, sigma, log=FALSE),
          (1/sigma) * dt((x-mu)/sigma, df=nu, log=FALSE))
     else dens <- if(nu > 1000000) dnorm(x, mu, sigma, log=FALSE)
     else (1/sigma) * dt((x-mu)/sigma, df=nu, log=FALSE)
     dens <- if(log == FALSE) dens else log(dens)
     return(dens)
     }
pst <- function(q, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     if(any(sigma <= 0)) stop("sigma must be positive in pst().")
     if(any(nu <= 0)) stop("nu must be positive in pst().")
     if(length(nu)>1) cdf <- ifelse(nu>1000000,
          pnorm(q, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          pt((q-mu)/sigma, df=nu,   lower.tail=lower.tail, log.p=log.p))
     else cdf <- if(nu > 1000000) pnorm(q, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)
     else pt((q-mu)/sigma, df=nu, lower.tail=lower.tail, log.p=log.p)
     return(cdf)
     }
qst <- function(p, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     if(any(sigma <= 0)) stop("sigma must be positive in qst().")
     if(any(nu <= 0)) stop("nu must be positive in qst().")
     if(any(p < 0)|any(p > 1)) stop("p must be between 0 and 1 in qst().")
     if(length(nu)>1) q <- ifelse(nu > 1000000,
          qnorm(p, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          mu + sigma * qt(p, df=nu, lower.tail=lower.tail))
     else q <- if(nu > 1000000) qnorm(p, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)
     else mu + sigma * qt(p, df=nu, lower.tail=lower.tail)
     return(q)
     }
rst <- function(n, mu=0, sigma=1, nu=10)
     {
     if(any(sigma <= 0)) stop("sigma must be positive in rst().")
     if(any(nu <= 0)) stop("nu must be positive in rst().")
     if(any(n <= 0)) stop("n must be a positive integer in rst().")
     n <- ceiling(n)
     p <- runif(n)
     r <- qst(p, mu=mu, sigma=sigma, nu=nu)
     return(r)
     }

###########################################################################
# Truncated Distribution                                                  #
#                                                                         #
# These functions are similar to those from Nadarajah, S. and Kotz, S.    #
# (2006). R Programs for Computing Truncated Distributions. Journal of    #
# Statistical Software, 16, Code Snippet 2, 1-8.                          #
###########################################################################

dtrunc <- function(x, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in dtrunc().")
     tt <- rep(0, length(x))
     g <- get(paste("d", spec, sep = ""), mode = "function")
     G <- get(paste("p", spec, sep = ""), mode = "function")
     tt[x>=a & x<=b] <- g(x[x>=a&x<=b], ...) / (G(b, ...) - G(a, ...))
     return(tt)
     }

extrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in extrunc().")
     f <- function(x) x * dtrunc(x, spec, a=a, b=b, ...)
     return(integrate(f, lower=a, upper=b)$value)
     }

ptrunc <- function(x, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in ptrunc().")
     tt <- x
     aa <- rep(a, length(x))
     bb <- rep(b, length(x))
     G <- get(paste("p", spec, sep = ""), mode = "function")
     tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
     tt <- tt - G(aa, ...)
     tt <- tt/(G(bb, ...) - G(aa, ...))
     return(tt)
     }

qtrunc <- function(p, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in qtrunc().")
     tt <- p
     G <- get(paste("p", spec, sep = ""), mode = "function")
     Gin <- get(paste("q", spec, sep = ""), mode = "function")
     tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
     return(tt)
     }

rtrunc <- function(n, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in rtrunc().")
     x <- u <- runif(n, min = 0, max = 1)
     x <- qtrunc(u, spec, a = a, b = b,...)
     return(x)
     }

vartrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b in vartrunc().")
     ex <- extrunc(spec, a = a, b = b, ...)
     f <- function(x) (x - ex)^2 * dtrunc(x, spec, a = a, b = b, ...)
     tt <- integrate(f, lower = a, upper = b)$value
     return(tt)
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
