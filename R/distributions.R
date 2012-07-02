###########################################################################
# Asymmetric Laplace Distribution                                         #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dalaplace <- function(x, location=0, scale=1, kappa=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(location), length(scale), length(kappa))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     logconst <- 0.5 * log(2) - log(scale) + log(kappa) - log1p(kappa^2)
     exponent <- -(sqrt(2) / scale) * abs(x - location) *
          ifelse(x >= location, kappa, 1/kappa)
     dens <- logconst + exponent
     if(log == FALSE) dens <- exp(logconst + exponent)
     return(dens)
     }
palaplace <- function(q, location=0, scale=1, kappa=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if((kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(location), length(scale), length(kappa))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     exponent <- -(sqrt(2) / scale) * abs(q - location) *
          ifelse(q >= location, kappa, 1/kappa)
     temp <- exp(exponent) / (1 + kappa^2)
     p <- 1 - temp
     index1 <- (q < location)
     p[index1] <- (kappa[index1])^2 * temp[index1]
     return(p)
     }
qalaplace <- function(p, location=0, scale=1, kappa=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(location), length(scale), length(kappa))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     q <- p
     temp <- kappa^2 / (1 + kappa^2)
     index1 <- (p <= temp)
     exponent <- p[index1] / temp[index1]
     q[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
          log(exponent) / sqrt(2)
     q[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
          (log1p((kappa[!index1])^2) + log1p(-p[!index1])) / sqrt(2)
     q[p == 0] = -Inf
     q[p == 1] = Inf
     return(q)
     }
ralaplace <- function(n, location=0, scale=1, kappa=1)
     {
     location <- rep(location, len=n)
     scale <- rep(scale, len=n)
     kappa <- rep(kappa, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     x <- location + scale *
          log(runif(n)^kappa / runif(n)^(1/kappa)) / sqrt(2)
     return(x)
     }

###########################################################################
# Asymmetric Log-Laplace Distribution                                     #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dallaplace <- function(x, location=0, scale=1, kappa=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(location), length(scale), length(kappa))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     exponent <- ifelse(x >= Delta, -(Alpha+1),
          (Beta-1)) * (log(x) - location)
     dens <- -location + log(Alpha) + log(Beta) -
          log(Alpha + Beta) + exponent
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pallaplace <- function(q, location=0, scale=1, kappa=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(location), length(scale), length(kappa))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     temp <- Alpha + Beta
     p <- (Alpha / temp) * (q / Delta)^(Beta)
     p[q <= 0] <- 0
     index1 <- (q >= Delta)
     p[index1] <- (1 - (Beta/temp) * (Delta/q)^(Alpha))[index1]
     return(p)
     }
qallaplace <- function(p, location=0, scale=1, kappa=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale); kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(location), length(scale), length(kappa))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN); kappa <- rep(kappa, len=NN)
     Alpha <- sqrt(2) * kappa / scale
     Beta  <- sqrt(2) / (scale * kappa)
     Delta <- exp(location)
     temp <- Alpha + Beta
     q <- Delta * (p * temp / Alpha)^(1/Beta)
     index1 <- (p > Alpha / temp)
     q[index1] <- (Delta * ((1-p) * temp / Beta)^(-1/Alpha))[index1]
     q[p == 0] <- 0
     q[p == 1] <- Inf
     return(q)
     }
rallaplace <- function(n, location=0, scale=1, kappa=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     kappa <- rep(kappa, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     x <- exp(location) *
          (runif(n)^kappa / runif(n)^(1/kappa))^(scale / sqrt(2))
     return(x)
     }

###########################################################################
# Bernoulli Distribution                                                  #
#                                                                         #
# These functions are similar to those in the Rlab package.               #
###########################################################################

dbern <- function(x, prob, log=FALSE)
     {return(dbinom(x, 1, prob, log))}
pbern <- function(q, prob, lower.tail=TRUE, log.p=FALSE)
     {return(pbinom(q, 1, prob, lower.tail, log.p))}
qbern <- function(p, prob, lower.tail=TRUE, log.p=FALSE)
     {return(qbinom(p, 1, prob, lower.tail, log.p))}
rbern <- function(n, prob)
     {return(rbinom(n, 1, prob))}

###########################################################################
# Categorical Distribution                                                #
###########################################################################

dcat <- function(x, p, log=FALSE)
     {
     if(is.vector(x) & !is.matrix(p))
          p <- matrix(p, length(x), length(p), byrow=TRUE)
     if(is.matrix(x) & !is.matrix(p))
          p <- matrix(p, nrow(x), length(p), byrow=TRUE)
     if(is.vector(x) & {length(x) == 1}) {
          temp <- rep(0, ncol(p))
          temp[x] <- 1
          x <- t(temp)}
     else if(is.vector(x) & (length(x) > 1))
          x <- as.indicator.matrix(x)
     if(!identical(nrow(x), nrow(p)))
          stop("The number of rows of x and p differ.")
     if(!identical(ncol(x), ncol(p))) {
          x.temp <- matrix(0, nrow(p), ncol(p))
          x.temp[,as.numeric(colnames(x))] <- x
          x <- x.temp}
     dens <- x*p
     if(log == TRUE) dens <- x*log(p)
     dens <- as.vector(rowSums(dens))
     return(dens)
     }
rcat <- function(n, p)
     {
     x <- which(rmultinom(n, size=1, prob=p) == 1, arr.ind=TRUE)[,'row']
     return(x)
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
          s <- sum({alpha-1}*log(x))
          exp(sum(s)-logD)
          }
    if(!is.matrix(x)) x <- matrix(x)
    if(is.data.frame(x)) x <- matrix(x)
    else x <- t(x)
    if(!is.matrix(alpha))
         alpha <- matrix(alpha, ncol=length(alpha), nrow=nrow(x),
              byrow=TRUE)
    if(!identical(dim(x), dim(alpha)))
         stop("The dimensions of x and alpha differ.")
    dens <- vector(length=nrow(x))
    for (i in 1:nrow(x)) {dens[i] <- dirichlet1(x[i,], alpha[i,])}
    dens[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    dens[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
    if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
    return(dens)
    }
rdirichlet <- function(n, alpha)
     {
     l <- length(alpha)
     x <- matrix(rgamma(l*n,alpha), ncol=l, byrow=TRUE)
     sm <- x %*% rep(1,l)
     return(x / as.vector(sm))
     }

###########################################################################
# Half-Cauchy Distribution                                                #
###########################################################################

dhalfcauchy <- function(x, scale=25, log=FALSE)
     {
     x <- as.vector(x); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(scale))
     x <- rep(x, len=NN); scale <- rep(scale, len=NN)
     dens <- 2*scale / (pi*{x*x + scale*scale})
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
phalfcauchy <- function(q, scale=25)
     {
     q <- as.vector(q); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     z <- {2/pi}*atan(q/scale)
     return(z)
     }
qhalfcauchy <- function(p, scale=25)
     {
     p <- as.vector(p); scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     q <- scale*tan({pi*p}/2)
     return(q)
     }
rhalfcauchy <- function(n, scale=25)
     {
     scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     p <- runif(n, 0, 1)
     x <- scale*tan({pi*p}/2)
     return(x)
     }

###########################################################################
# Half-Normal Distribution                                                #
#                                                                         #
# This half-normal distribution has mean=0 and is similar to the halfnorm #
# functions in package fdrtool.                                           #
###########################################################################

dhalfnorm <- function(x, scale=sqrt(pi/2), log=FALSE)
     {
     dens <- 2*dnorm(x, mean=0, sd=sqrt(pi/2) / scale)
     if(log == TRUE) dens <- log(dens  + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
phalfnorm <- function(q, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     p <- 2*pnorm(q, mean=0, sd=sqrt(pi/2) / scale) - 1
     if(lower.tail == FALSE) p <- 1-p
     if(log.p == TRUE) p <- log.p(p)
     return(p)
     }
qhalfnorm <- function(p, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     if(log.p == TRUE) p <- exp(p)
     if(lower.tail == FALSE) p <- 1-p
     q <- qnorm((p+1)/2, mean=0, sd=sqrt(pi/2) / scale)
     return(q)
     }
rhalfnorm <- function(n, scale=sqrt(pi/2))
     {
     scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- abs(rnorm(n, mean=0, sd=sqrt(pi/2) / scale))
     return(x)
     }

###########################################################################
# Half-t Distribution                                                     #
###########################################################################

dhalft <- function(x, scale=25, nu=1, log=FALSE)
     {
     x <- as.vector(x); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(scale), length(nu))
     x <- rep(x, len=NN); scale <- rep(scale, len=NN)
     nu <- rep(nu, len=NN)
     dens <- (1 + {1/nu}*{x/scale}*{x/scale})^(-{nu+1}/2)
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
phalft <- function(q, scale=25, nu=1)
     {
     q <- as.vector(q); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(scale), length(nu))
     q <- rep(q, len=NN); scale <- rep(scale, len=NN)
     p <- ptrunc(q, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(p)
     }
qhalft <- function(p, scale=25, nu=1)
     {
     p <- as.vector(p); scale <- as.vector(scale); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(scale), length(nu))
     p <- rep(p, len=NN); scale <- rep(scale, len=NN)
     q <- rtrunc(p, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(q)
     }
rhalft <- function(n, scale=25, nu=1)
     {
     scale <- rep(scale, len=n); nu <- rep(nu, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- rtrunc(n, "st", a=0, b=Inf, mu=0, sigma=scale, nu=nu)
     return(x)
     }

###########################################################################
# Inverse Chi-Squared Distribution                                        #
#                                                                         #
# These functions are similar to those in the GeoR package.               #
###########################################################################

dinvchisq <- function(x, df, scale=1/df, log=FALSE)
     {
     x <- as.vector(x); df <- as.vector(df); scale <- as.vector(scale)
     if(any(x <= 0)) stop("x must be positive.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(df), length(scale))
     x <- rep(x, len=NN); df <- rep(df, len=NN);
     scale <- rep(scale, len=NN)
     nu <- df / 2
     dens <- nu*log(nu) - log(gamma(nu)) + nu*log(scale) -
          (nu+1)*log(x) - (nu*scale/x)
     if(log == FALSE) dens <- exp(dens)
     }

rinvchisq <- function(n, df, scale=1/df)
     {
     df <- rep(df, len=n); scale <- rep(scale, len=n)
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     z <- rchisq(n, df=df)
     z <- ifelse(z == 0, 1e-100, z)
     x <- (df*scale) / z
     return(x)
     }

###########################################################################
# Inverse Gamma Distribution                                              #
#                                                                         #
# These functions are similar to those in the MCMCpack package.           #
###########################################################################

dinvgamma <- function(x, shape=1, scale=1, log=FALSE)
     {
     x <- as.vector(x); shape <- as.vector(shape)
     scale <- as.vector(scale)
     if(any(shape <= 0) | any(scale <=0))
          stop("The shape and scale parameters must be positive.")
     NN <- max(length(x), length(shape), length(scale))
     x <- rep(x, len=NN); shape <- rep(shape, len=NN)
     scale <- rep(scale, len=NN)
     alpha <- shape; beta <- scale
     dens <- alpha * log(beta) - lgamma(alpha) -
          {alpha + 1} * log(x) - {beta/x}
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rinvgamma <- function(n, shape=1, scale=1)
     {return(1 / rgamma(n=n, shape=shape, rate=scale))}

###########################################################################
# Inverse Gaussian Distribution                                           #
###########################################################################

dinvgaussian <- function(x, mu, lambda, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); lambda <- as.vector(lambda)
     if(any(x <= 0)) stop("x must be positive.")
     if(any(mu <= 0)) stop("The mu parameter must be positive.")
     if(any(lambda <= 0)) stop("The lambda parameter must be positive.")
     NN <- max(length(x), length(mu), length(lambda))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     lambda <- rep(lambda, len=NN)
     dens <- (lambda / (2*pi*x^3))^0.5 *
          exp(-((lambda*(x - mu)^2) / (2*mu^2*x)))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rinvgaussian <- function(n, mu, lambda)
     {
     mu <- rep(mu, len=n); lambda <- rep(lambda, len=n)
     if(any(mu <= 0)) stop("The mu parameter must be positive.")
     if(any(lambda <= 0)) stop("The lambda parameter must be positive.")
     nu <- rnorm(n)
     y <- nu^2
     x <- mu + ((mu^2*y)/(2*lambda)) - (mu/(2*lambda)) *
          sqrt(4*mu*lambda*y + mu^2*y^2)
     z <- runif(n)
     x <- ifelse(z > (mu / (mu+x)), mu^2/x, x)
     return(x)
     }

###########################################################################
# Inverse Wishart Distribution                                            #
###########################################################################

dinvwishart <- function(Sigma, nu, S, log=FALSE)
     {
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(!identical(dim(S), dim(Sigma)))
          stop("The dimensions of Sigma and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Sigma)
     detS <- det(S)
     detSigma <- det(Sigma)
     invSigma <- as.inverse(Sigma)
     gamprod <- 1
     for (i in 1:k) {gamprod <- gamprod * gamma((nu + 1 - i) / 2)}
     dens <- (2^(nu*k/2)*pi^(k*(k-1)/4)*gamprod)^(-1) * detS^(nu/2) *
          detSigma^(-(nu+k+1)/2) * exp(-0.5 * tr(S %*% invSigma))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rinvwishart <- function(nu, S)
     {return(as.inverse(rwishart(nu, as.inverse(S))))}

###########################################################################
# Laplace Distribution                                                    #
#                                                                         #
# These functions are similar to those in the VGAM package.               #
###########################################################################

dlaplace <- function(x, location=0, scale=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(location), length(scale))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     dens <- (-abs(x - location) / scale) - log(2 * scale)
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
plaplace <- function(q, location=0, scale=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     z <- {q - location} / scale
     NN <- max(length(q), length(location), length(scale))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     p <- ifelse(q < location, 0.5 * exp(z), 1 - 0.5 * exp(-z))
     return(p)
     }
qlaplace <- function(p, location=0, scale=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(location), length(scale))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     q <- location - sign(p - 0.5) * scale * log(2 * ifelse(p < 0.5,
          p, 1 - p))
     return(q)
     }
rlaplace <- function(n, location=0, scale=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     r <- runif(n)
     x <- location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5,
          r, 1 - r))
     return(x)
     }

###########################################################################
# Laplace Distribution (Precision Parameterization)                       #
###########################################################################

dlaplacep <- function(x, mu=0, tau=1, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     dens <- (tau/2) * exp(-tau*abs(x-mu))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
plaplacep <- function(q, mu=0, tau=1)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     p <- plaplace(q, mu, scale=1/tau)
     return(p)
     }
qlaplacep <- function(p, mu=0, tau=1)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     q <- qlaplace(p, mu, scale=1/tau)
     return(q)
     }
rlaplacep <- function(n, mu=0, tau=1)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     x <- rlaplace(n, mu, scale=1/tau)
     return(x)
     }

###########################################################################
# Log-Laplace Distribution                                                #
###########################################################################

dllaplace <- function(x, location=0, scale=1, log=FALSE)
     {
     x <- as.vector(x); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(x), length(location), length(scale))
     x <- rep(x, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     exponent <- ifelse(x >= Delta, -(Alpha+1),
          (Beta-1)) * (log(x) - location)
     dens <- -location + log(Alpha) + log(Beta) -
          log(Alpha + Beta) + exponent
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
pllaplace <- function(q, location=0, scale=1)
     {
     q <- as.vector(q); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(q), length(location), length(scale))
     q <- rep(q, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     temp <- Alpha + Beta
     p <- (Alpha / temp) * (q / Delta)^(Beta)
     p[q <= 0] <- 0
     index1 <- (q >= Delta)
     p[index1] <- (1 - (Beta/temp) * (Delta/q)^(Alpha))[index1]
     return(p)
     }
qllaplace <- function(p, location=0, scale=1)
     {
     p <- as.vector(p); location <- as.vector(location)
     scale <- as.vector(scale)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     NN <- max(length(p), length(location), length(scale))
     p <- rep(p, len=NN); location <- rep(location, len=NN)
     scale <- rep(scale, len=NN)
     Alpha <- sqrt(2) * scale
     Beta  <- sqrt(2) / scale
     Delta <- exp(location)
     temp <- Alpha + Beta
     q <- Delta * (p * temp / Alpha)^(1/Beta)
     index1 <- (p > Alpha / temp)
     q[index1] <- (Delta * ((1-p) * temp / Beta)^(-1/Alpha))[index1]
     q[p == 0] <- 0
     q[p == 1] <- Inf
     return(q)
     }
rllaplace <- function(n, location=0, scale=1)
     {
     location <- rep(location, len=n); scale <- rep(scale, len=n)
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     x <- exp(location) * (runif(n) / runif(n))^(scale / sqrt(2))
     return(x)
     }

###########################################################################
# Log-Normal Distribution (Precision Parameterization)                    #
###########################################################################

dlnormp <- function(x, mu, tau, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     dens <- sqrt(tau/(2*pi)) * (1/x) * exp(-(tau/2)*(log(x-mu))^2)
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
plnormp <- function(q, mu, tau, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     p <- pnorm(q, mu, sqrt(1/tau), lower.tail, log.p)
     return(p)
     }
qlnormp <- function(p, mu, tau, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu); tau <- as.vector(tau)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN); tau <- rep(tau, len=NN)
     q <- qnorm(p, mu, sqrt(1/tau), lower.tail, log.p)
     return(q)
     }
rlnormp <- function(n, mu, tau)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     x <- rnorm(n, mu, sqrt(1/tau))
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution                                        #
###########################################################################

dmvc <- function(x, mu, S, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(S)) S <- diag(ncol(x))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     k <- nrow(S)
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(gamma(k/2) / (gamma(0.5) * 1^(k/2) *
          pi^(k/2) * sqrt(det(S)) * (1 + z)^((1+k)/2)))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rmvc <- function(n=1, mu=rep(0,k), S)
     {
     mu <- as.vector(mu)
     if(missing(S)) S <- diag(length(mu))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     k <- ncol(S)
     x <- rchisq(n,1)
     x <- ifelse(x == 0, 1e-100, x)
     z <- rmvn(n, rep(0,k), S)
     x <- t(mu + t(z/sqrt(x)))
     return(x)
     }

###########################################################################
# Multivariate Cauchy Distribution (Precision Parameterization)           #
###########################################################################

dmvcp <- function(x, mu, Omega, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((gamma((1+k)/2) /
          (gamma(0.5)*1^(k/2)*pi^(k/2))) * detOmega^0.5 *
          (1 + z)^(-(1+k)/2))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rmvcp <- function(n=1, mu, Omega)
     {
     mu <- as.vector(mu)
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     Sigma <- as.inverse(Omega)
     k <- ncol(Sigma)
     x <- rchisq(n,1)
     x <- ifelse(x == 0, 1e-100, x)
     z <- rmvn(n, rep(0,k), Sigma)
     x <- t(mu + t(z/sqrt(x)))
     return(x)
     }

###########################################################################
# Multivariate Laplace Distribution                                       #
###########################################################################

dmvl <- function(x, mu, Sigma, log=FALSE)
     {
     return(dmvpe(x, mu, Sigma, kappa=0.5, log))
     }
rmvl <- function(n, mu, Sigma)
     {
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     mu <- matrix(mu, n, ncol(Sigma), byrow=TRUE)
     e <- matrix(rexp(n, 1), n, ncol(mu))
     z <- rmvn(n, rep(0, ncol(mu)), Sigma)
     x <- e*mu + sqrt(e)*z
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution                                        #
###########################################################################

dmvn <- function(x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- nrow(Sigma)
     Omega <- as.inverse(Sigma)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     d <- eigen(Sigma, symmetric = TRUE)$values
     dens <- as.vector(-0.5 * (k * log(2 * pi) + sum(log(d))) - (0.5 * z))
     if(log == FALSE) dens <- exp(dens)
     return(dens)
     }
rmvn <- function(n=1, mu=rep(0,k), Sigma)
     {
     mu <- as.vector(mu)
     if(missing(Sigma)) Sigma <- diag(length(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     k <- ncol(Sigma)
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     x <- t(mu + t(z))
     return(x)
     }

###########################################################################
# Multivariate Normal Distribution (Precision Parameterization)           #
###########################################################################

dmvnp <- function(x, mu, Omega, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     k <- nrow(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((2*pi)^(-k/2) * detOmega^0.5 * exp(-0.5 * z))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rmvnp <- function(n=1, mu=rep(0,k), Omega)
     {
     mu <- as.vector(mu)
     if(missing(Omega)) Omega <- diag(length(mu))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     Sigma <- as.inverse(Omega)
     k <- ncol(Sigma)
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     x <- t(mu + t(z))
     return(x)
     }

###########################################################################
# Multivariate Polya Distribution                                         #
###########################################################################

dmvpolya <- function(x, alpha, log=FALSE)
     {
     x <- as.vector(x)
     alpha <- as.vector(alpha)
     if(!identical(length(x), length(alpha)))
          stop("x and alpha differ in length.")
     dens <- (factorial(sum(x)) / prod(factorial(x))) *
          (factorial(sum(alpha)-1) / factorial(sum(x) + sum(alpha)-1)) *
          (prod(factorial(x + alpha - 1) / factorial(alpha - 1)))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }

###########################################################################
# Multivariate Power Exponential Distribution                             #
###########################################################################

dmvpe <- function(x=c(0,0), mu=c(0,0), Sigma=diag(2), kappa=1, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     k <- nrow(Sigma)
     ss <- x - mu
     temp <- rowSums({ss %*% Sigma} * ss)
     dens <- ((k*gamma(k/2)) / (pi^(k/2) * sqrt(det(Sigma)) *
          gamma(1 + k/(2*kappa)) * 2^(1 + k/(2*kappa)))) *
          exp(-0.5*temp)^kappa
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(as.vector(dens))
     }

###########################################################################
# Multivariate t Distribution                                             #
###########################################################################

dmvt <- function(x, mu, S, df=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(S)) S <- diag(ncol(x))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     if(any(df > 10000)) return(dmvn(x, mu, S, log))
     k <- nrow(S)
     ss <- x - mu
     Omega <- as.inverse(S)
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector(gamma((df+k)/2) / (gamma(df/2) * df^(k/2) *
          pi^(k/2) * sqrt(det(S)) * (1 + (1/df) * z)^((df+k)/2)))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rmvt <- function(n=1, mu=rep(0,k), S, df=Inf)
     {
     mu <- as.vector(mu)
     if(missing(S)) S <- diag(length(mu))
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(any(df <= 0)) stop("The df parameter must be positive.")
     k <- ncol(S)
     if(df==Inf) x <- 1 else x <- rchisq(n,df) / df
     x <- ifelse(x == 0, 1e-100, x)
     z <- rmvn(n, rep(0,k), S)
     x <- t(mu + t(z/sqrt(x)))
     return(x)
     }

###########################################################################
# Multivariate t Distribution (Precision Parameterization)                #
###########################################################################

dmvtp <- function(x, mu, Omega, nu=Inf, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rep(mu, each=nrow(x))
     if(missing(Omega)) Omega <- diag(ncol(x))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     if(any(nu > 10000)) return(dmvnp(x, mu, Omega, log))
     k <- ncol(Omega)
     detOmega <- det(Omega)
     ss <- x - mu
     z <- rowSums({ss %*% Omega} * ss)
     dens <- as.vector((gamma((nu+k)/2) /
          (gamma(nu/2)*nu^(k/2)*pi^(k/2))) * detOmega^0.5 *
          (1 + (1/nu) * z)^(-(nu+k)/2))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rmvtp <- function(n=1, mu, Omega, nu=Inf)
     {
     mu <- as.vector(mu)
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     Sigma <- as.inverse(Omega)
     k <- ncol(Sigma)
     if(nu == Inf) x <- 1 else x <- rchisq(n,nu) / nu
     x <- ifelse(x == 0, 1e-100, x)
     z <- rmvn(n, rep(0,k), Sigma)
     x <- t(mu + t(z/sqrt(x)))
     return(x)
     }

###########################################################################
# Normal Distribution (Precision Parameterization)                        #
###########################################################################

dnormp <- function(x, mean=0, prec=1, log=FALSE)
     {
     #dens <- sqrt(prec/(2*pi)) * exp(-(prec/2)*(x-mu)^2)
     dens <- dnorm(x, mean, sqrt(1/prec), log)
     return(dens)
     }
pnormp <- function(q, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
     {return(pnorm(q, mean=mean, sd=sqrt(1/prec), lower.tail, log.p))}
qnormp <- function(p, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
     {return(qnorm(p, mean=mean, sd=sqrt(1/prec), lower.tail, log.p))}
rnormp <- function(n, mean=0, prec=1)
     {return(rnorm(n, mean=mean, sd=sqrt(1/prec)))}

###########################################################################
# Normal Distribution (Variance Parameterization)                         #
###########################################################################

dnormv <- function(x, mean=0, var=1, log=FALSE)
     {
     #dens <- (1/(sqrt(2*pi*var))) * exp(-((x-mu)^2/(2*var)))
     dens <- dnorm(x, mean, sqrt(var), log)
     return(dens)
     }
pnormv <- function(q, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
     {return(pnorm(q, mean=mean, sd=sqrt(var), lower.tail, log.p))}
qnormv <- function(p, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
     {return(qnorm(p, mean=mean, sd=sqrt(var), lower.tail, log.p))}
rnormv <- function(n, mean=0, var=1)
     {return(rnorm(n, mean=mean, sd=sqrt(var)))}

###########################################################################
# Pareto Distribution                                                     #
###########################################################################

dpareto <- function(x, alpha, log=FALSE)
     {
     x <- as.vector(x); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(x), length(alpha))
     x <- rep(x, len=NN); alpha <- rep(alpha, len=NN)
     dens <- ifelse(x < 1, 0, alpha / x^(alpha + 1))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }

ppareto <- function(q, alpha)
     {
     q <- as.vector(q); alpha <- as.vector(alpha)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(q), length(alpha))
     q <- rep(q, len=NN); alpha <- rep(alpha, len=NN)
     p <- ifelse(q < 1, 0, 1 - 1/q^alpha)
     return(p)
     }

qpareto <- function(p, alpha)
     {
     p <- as.vector(p); alpha <- as.vector(alpha)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     NN <- max(length(p), length(alpha))
     p <- rep(p, len=NN); alpha <- rep(alpha, len=NN)
     q <- (1-p)^(-1/alpha)
     return(q)
     }

rpareto <- function(n, alpha)
     {
     alpha <- rep(alpha, len=n)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     x <- runif(n)^(-1/alpha)
     return(x)
     }

###########################################################################
# Power Exponential Distribution                                          #
#                                                                         #
# These functions are similar to those in the normalp package.            #
###########################################################################

dpe <- function(x, mu=0, sigma=1, kappa=2, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(x), length(mu), length(sigma), length(kappa))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     cost <- 2 * kappa^(1/kappa) * gamma(1 + 1/kappa) * sigma
     expon1 <- (abs(x - mu))^kappa
     expon2 <- kappa * sigma^kappa
     dens <- (1/cost) * exp(-expon1 / expon2)
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
ppe <- function(q, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(q), length(mu), length(sigma), length(kappa))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     z <- (q - mu) / sigma
     zz <- abs(z)^kappa
     p <- pgamma(zz, shape=1/kappa, scale=kappa)
     p <- p / 2
     p <- ifelse(z < 0, 0.5 - p, 0.5 + p)
     if(lower.tail == FALSE) p <- 1 - p
     if(log.p == TRUE) p <- log(p)
     return(p)
     }
qpe <- function(p, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu); sigma <- as.vector(sigma)
     kappa <- as.vector(kappa)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     NN <- max(length(p), length(mu), length(sigma), length(kappa))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); kappa <- rep(kappa, len=NN)
     if(log.p == TRUE) p <- log(p)
     if(lower.tail == FALSE) p <- 1 - p
     zp <- ifelse(p < 0.5, 0.5 - p, p - 0.5)
     zp <- 2 * zp
     qg <- qgamma(zp, shape=1/kappa, scale=kappa)
     z <- qg^(1/kappa)
     z <- ifelse(p < 0.5, -z, z)
     q <- mu + z * sigma
     return(q)
     }
rpe <- function(n, mu=0, sigma=1, kappa=2)
     {
     mu <- rep(mu, len=n); sigma <- rep(sigma, len=n)
     kappa <- rep(kappa, len=n)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(kappa <= 0)) stop("The kappa parameter must be positive.")
     qg <- rgamma(n, shape=1/kappa, scale=kappa)
     z <- qg^(1/kappa)
     z <- ifelse(runif(n) < 0.5, -z, z)
     x <- mu + z * sigma
     return(x)
     }

###########################################################################
# Skew Discrete Laplace Distribution                                      #
###########################################################################

dsdlaplace <- function(x, p, q, log=FALSE)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(x), length(p), length(q))
     x <- rep(x, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     dens <- ifelse(x >= 0, (1-p)*(1-q)/(1-p*q)*p^x,
          (1-p)*(1-q)/(1-p*q)*q^abs(x))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
psdlaplace <- function(x, p, q)
     {
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(x), length(p), length(q))
     x <- rep(x, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     pr <- ifelse(x >= 0, 1-(1-q)*p^(floor(x)+1)/(1-p*q),
          (1-p)*q^(-floor(x))/(1-p*q))
     return(pr)
     }
qsdlaplace <- function(prob, p, q)
     {
     if(any(prob < 0) || any(prob > 1)) stop("prob must be in [0,1].")
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(q < 0) || any(q > 1)) stop("q must be in [0,1].")
     NN <- max(length(prob), length(p), length(q))
     prob <- rep(prob, len=NN); p <- rep(p, len=NN); q <- rep(q, len=NN)
     x <- numeric(NN)
     for (i in 1:NN) {
          k <- 0
          if(prob[i] >= psdlaplace(k, p[i], q[i])) {
               while(prob[i] >= psdlaplace(k, p[i], q[i])) {
                    k <- k + 1}}
          else if(prob[i] < psdlaplace(k, p[i], q[i])) {
               while(prob[i] < psdlaplace(k, p[i], q[i])) {
                    k <- k - 1}
               k <- k + 1}
          x[i] <- k
          }
     return(x)
     }
rsdlaplace <- function(n, p, q)
     {
     if(length(p) > 1) stop("p must have a length of 1.")
     if(length(q) > 1) stop("q must have a length of 1.")
     if((p < 0) || (p > 1)) stop("p must be in [0,1].")
     if((q < 0) || (q > 1)) stop("q must be in [0,1].")
     u <- runif(n)
     return(qsdlaplace(u,p,q))
     }

###########################################################################
# Skew-Laplace Distribution                                               #
###########################################################################

dslaplace <- function(x, mu, alpha, beta, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(x), length(mu), length(alpha), length(beta))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     belowMu <- (1/ab) * exp((x - mu)/alpha)
     aboveMu <- (1/ab) * exp((mu - x)/beta)
     dens <- ifelse(x <= mu, belowMu, aboveMu)
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
pslaplace <- function(q, mu, alpha, beta)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(q), length(mu), length(alpha), length(beta))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     belowMu <- (alpha/ab) * exp((q - mu)/alpha)
     aboveMu <- 1 - (beta/ab) * exp((mu - q)/beta)
     p <- ifelse(q < mu, belowMu, aboveMu)
     return(p)
     }
qslaplace <- function(p, mu, alpha, beta)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     alpha <- as.vector(alpha); beta <- as.vector(beta)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     NN <- max(length(p), length(mu), length(alpha), length(beta))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     alpha <- rep(alpha, len=NN); beta <- rep(beta, len=NN)
     ab <- alpha + beta
     belowMu <- alpha*log(p*ab/alpha) + mu
     aboveMu <- mu - beta*log(ab*(1 - p)/beta)
     q <- ifelse(p < alpha/ab, belowMu, aboveMu)
     return(q)
     }
rslaplace <- function(n, mu, alpha, beta)
     {
     #mu <- rep(mu, len=n); alpha <- rep(alpha, len=n)
     #beta <- rep(beta, len=n)
     if(any(alpha <= 0)) stop("The alpha parameter must be positive.")
     if(any(beta <= 0)) stop("The beta parameter must be positive.")
     ab <- alpha + beta
     y <- rexp(n,1)
     probs <- c(alpha,beta) / ab
     signs <- sample(c(-1,1), n, replace=TRUE, prob=probs)
     mult <- ifelse(signs < 0, signs*alpha, signs*beta)
     x <- mult*y + mu
     return(x)
     }

###########################################################################
# Stick-Breaking Prior Distribution                                       #
###########################################################################

dStick <- function(theta, gamma, log=FALSE)
  {
  dens <- sum(dbeta(theta, 1, gamma, log=log))
  return(dens)
  }

###########################################################################
# Student t Distribution (3-parameter)                                    #
#                                                                         #
# These functions are similar to those in the gamlss.dist package.        #
###########################################################################

dst <- function(x, mu=0, sigma=1, nu=10, log=FALSE)
     {return({1/sigma} * dt({x-mu}/sigma, df=nu, log))}
pst <- function(q, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     sigma <- as.vector(sigma); nu <- as.vector(nu)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(q), length(mu), length(sigma), length(nu))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); nu <- rep(nu, len=NN)
     if(length(nu) > 1) cdf <- ifelse(nu > 1000000,
          pnorm(q, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          pt({q-mu}/sigma, df=nu, lower.tail=lower.tail, log.p=log.p))
     else cdf <- if(nu > 1000000) {pnorm(q, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)}
          else {pt({q-mu}/sigma, df=nu, lower.tail=lower.tail,
               log.p=log.p)}
     return(cdf)
     }
qst <- function(p, mu=0, sigma=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     sigma <- as.vector(sigma); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(p), length(mu), length(sigma), length(nu))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     sigma <- rep(sigma, len=NN); nu <- rep(nu, len=NN)
     if(length(nu) > 1) q <- ifelse(nu > 1000000,
          qnorm(p, mu, sigma, lower.tail=lower.tail, log.p=log.p),
          mu + sigma * qt(p, df=nu, lower.tail=lower.tail))
     else q <- if(nu > 1000000) {qnorm(p, mu, sigma,
          lower.tail=lower.tail, log.p=log.p)}
          else {mu + sigma * qt(p, df=nu, lower.tail=lower.tail)}
     return(q)
     }
rst <- function(n, mu=0, sigma=1, nu=10)
     {
     mu <- rep(mu, len=n); sigma <- rep(sigma, len=n); nu <- rep(nu, len=n)
     if(any(sigma <= 0)) stop("The sigma parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     n <- ceiling(n)
     p <- runif(n)
     x <- qst(p, mu=mu, sigma=sigma, nu=nu)
     return(x)
     }

###########################################################################
# Student t Distribution (Precision Parameterization)                     #
###########################################################################

dstp <- function(x, mu=0, tau=1, nu=10, log=FALSE)
     {
     x <- as.vector(x); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(x), length(mu), length(tau), length(nu))
     x <- rep(x, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     if(length(nu > 1)) {
          dens <- ifelse(nu > 1000000,
               dnorm(x, mu, sqrt(1/tau), log=FALSE),
               (gamma((nu+1)/2) / gamma(nu/2)) * sqrt(tau/(nu*pi)) *
                    (1 + (tau/nu)*(x-mu)^2)^(-(nu+1)/2))}
     else {dens <- if(nu > 1000000) {dnorm(x, mu, sqrt(1/tau), log=FALSE)}
          else {(gamma((nu+1)/2) / gamma(nu/2)) * sqrt(tau/(nu*pi)) *
                    (1 + (tau/nu)*(x-mu)^2)^(-(nu+1)/2)}}
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
pstp <- function(q, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     q <- as.vector(q); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(q), length(mu), length(tau), length(nu))
     q <- rep(q, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     cdf <- pst(q, mu, sqrt(1/tau), nu, lower.tail, log.p)
     return(cdf)
     }
qstp <- function(p, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
     {
     p <- as.vector(p); mu <- as.vector(mu)
     tau <- as.vector(tau); nu <- as.vector(nu)
     if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     NN <- max(length(p), length(mu), length(tau), length(nu))
     p <- rep(p, len=NN); mu <- rep(mu, len=NN)
     tau <- rep(tau, len=NN); nu <- rep(nu, len=NN)
     q <- qst(p, mu, sqrt(1/tau), nu, lower.tail, log.p)
     return(q)
     }
rstp <- function(n, mu=0, tau=1, nu=10)
     {
     mu <- rep(mu, len=n); tau <- rep(tau, len=n); nu <- rep(nu, len=n)
     if(any(tau <= 0)) stop("The tau parameter must be positive.")
     if(any(nu <= 0)) stop("The nu parameter must be positive.")
     x <- rst(n, mu, sqrt(1/tau), nu)
     n <- ceiling(n)
     p <- runif(n)
     x <- qst(p, mu=mu, sqrt(1/tau), nu=nu)
     return(x)
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
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     tt <- rep(0, length(x))
     g <- get(paste("d", spec, sep=""), mode="function")
     G <- get(paste("p", spec, sep=""), mode="function")
     tt[x>=a & x<=b] <- g(x[x>=a & x<=b], ...) / (G(b, ...) - G(a, ...))
     return(tt)
     }
extrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     f <- function(x) x * dtrunc(x, spec, a=a, b=b, ...)
     return(integrate(f, lower=a, upper=b)$value)
     }
ptrunc <- function(x, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     tt <- x
     aa <- rep(a, length(x))
     bb <- rep(b, length(x))
     G <- get(paste("p", spec, sep=""), mode="function")
     tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
     tt <- tt - G(aa, ...)
     tt <- tt / {G(bb, ...) - G(aa, ...)}
     return(tt)
     }
qtrunc <- function(p, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     tt <- p
     G <- get(paste("p", spec, sep=""), mode="function")
     Gin <- get(paste("q", spec, sep=""), mode="function")
     tt <- Gin(G(a, ...) + p*{G(b, ...) - G(a, ...)}, ...)
     return(tt)
     }
rtrunc <- function(n, spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     x <- u <- runif(n)
     x <- qtrunc(u, spec, a = a, b = b,...)
     return(x)
     }
vartrunc <- function(spec, a=-Inf, b=Inf, ...)
     {
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     ex <- extrunc(spec, a = a, b = b, ...)
     f <- function(x) {x - ex}^2 * dtrunc(x, spec, a = a, b = b, ...)
     tt <- integrate(f, lower = a, upper = b)$value
     return(tt)
     }

###########################################################################
# Wishart Distribution                                                    #
###########################################################################

dwishart <- function(Omega, nu, S, log=FALSE)
     {
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(!identical(dim(Omega), dim(S)))
          stop("The dimensions of Omega and S differ.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     k <- nrow(Omega)
     detS <- det(S)
     detOmega <- det(Omega)
     invS <- as.inverse(S)
     gamprod <- 1
     for (i in 1:k) {gamprod <- gamprod * gamma((nu + 1 - i) / 2)}
     dens <- (2^(nu*k/2)*pi^(k*(k-1)/4)*gamprod)^(-1) * detS^(-nu/2) *
          detOmega^((nu-k-1)/2) * exp(-0.5 * tr(invS %*% Omega))
     if(log == TRUE) dens <- log(dens + .Machine$double.xmin) #Prevent -Inf
     return(dens)
     }
rwishart <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.definite(S))
          stop("Matrix S is not positive-definite.")
     if(nu < nrow(S)) {
          stop("The nu parameter is less than the dimension of S.")}
     k <- nrow(S)
     Z <- matrix(0, k, k)
     x <- rchisq(k, nu:{nu - k + 1})
     x <- ifelse(x == 0, 1e-100, x)
     diag(Z) <- sqrt(x)
     if(k > 1) {
          kseq <- 1:(k-1)
          Z[rep(k*kseq, kseq) +
               unlist(lapply(kseq, seq))] <- rnorm(k*{k-1}/2)}
     return(crossprod(Z %*% chol(S)))
     }

#End
