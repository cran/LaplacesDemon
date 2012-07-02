###########################################################################
# LML (Logarithm of the Marginal Likelihood)                              #
#                                                                         #
# This function provides a few different methods of approximating the     #
# logarithm of the marginal likelihood (LML), and usually also            #
# approximates the Hessian matrix (returning a variance-covariance        #
# matrix, its negative inverse).                                          #
###########################################################################

LML <- function(Model=NULL, Data=NULL, Modes=NULL, theta=NULL, LL=NULL,
     method="NSIS")
     {
     LML.out <- list(LML=NA, VarCov=NA)
     if(!is.null(LL) & method == "HME") {
          med <- median(LL)
          LML <- med - log(mean(exp(-LL + med)))
          ### Output
          LML.out <- list(LML=LML, VarCov=NA)
          }
     if(!is.null(Model) & !is.null(Data) & !is.null(Modes) &
          (method == "LME1")) {
          Interval <- 1.0E-6
          parm.len <- length(Modes)
          eps <- Interval * Modes
          timetest <- system.time(Model(Modes, Data))
          timeest <- parm.len^2 * 2 * as.vector(timetest[3]) * 1.1 / 60
          cat("\nApproximating LML should take just over",
               round(timeest), "minutes.\n")
          Approx.Hessian <- Hessian(Model, Modes, Data)
          Inverse.test <- try(-as.inverse(Approx.Hessian), silent=TRUE)
          if(!inherits(Inverse.test, "try-error")) {
               VarCov <- Inverse.test
               diag(VarCov) <- ifelse(diag(VarCov) <= 0,
                    .Machine$double.eps, diag(VarCov))}
          else {
               cat("\nWARNING: Failure to solve matrix inversion of ",
                    "Approx. Hessian in LML.\n", sep="")
               cat("NOTE: Identity matrix is supplied instead.\n")
               VarCov <- diag(parm.len)}
          ### Logarithm of the Marginal Likelihood
          LML <- NA
          options(warn=-1)
          LML.test <- try(parm.len/2 * log(2*pi) + 0.5*log(det(VarCov)) +
               as.vector(Model(Modes, Data)[[1]]), silent=TRUE)
          if(!inherits(LML.test, "try-error")) LML <- LML.test[1]
          options(warn=0)
          ### Output
          LML.out <- list(LML=LML, VarCov=VarCov)
          }
     if(!is.null(Model) & !is.null(Data) & !is.null(Modes) &
          (method == "LME2")) {
          parm.len <- length(Modes)
          KoschalHess <- function(pars, fun, ...,
               .relStep=(.Machine$double.eps)^(1/3), minAbsPar=0)
               {
               pars <- as.numeric(pars)
               npar <- length(pars)
               incr <- ifelse(abs(pars) <= minAbsPar, minAbsPar * .relStep,
                    abs(pars) * .relStep)
               baseInd <- diag(npar)
               frac <- c(1, incr, incr^2)
               cols <- list(0, baseInd, -baseInd)
               for (i in seq_along(pars)[-npar]) {
                    cols <- c(cols, list(baseInd[, i] + baseInd[, -(1:i)]))
                    frac <- c(frac, incr[i] * incr[-(1:i)])
                    }
               indMat <- do.call("cbind", cols)
               shifted <- pars + incr * indMat
               indMat <- t(indMat)
               Xcols <- list(1, indMat, indMat^2)
               for (i in seq_along(pars)[-npar])
                    Xcols <- c(Xcols, list(indMat[, i] * indMat[, -(1:i)]))
               coefs <- solve(do.call("cbind", Xcols) ,
                    apply(shifted, 2, fun, ...) )/frac
               Hess <- diag(coefs[1 + npar + seq_along(pars)], ncol=npar)
               Hess[row(Hess) > col(Hess)] <- coefs[-(1:(1 + 2 * npar))]
               list(mean=coefs[1], gradient=coefs[1 + seq_along(pars)],
                    Hessian=(Hess + t(Hess)))
               }
          Koschal.test <- try(KoschalHess(Modes,
               function(x) Model(x, Data)[[1]])$Hessian, silent=TRUE)
          if(!inherits(Koschal.test, "try-error"))
               Approx.Hessian <- Koschal.test
          Inverse.test <- try(-as.inverse(Approx.Hessian), silent=TRUE)
          if(!inherits(Inverse.test, "try-error")) {
               VarCov <- Inverse.test
               diag(VarCov) <- ifelse(diag(VarCov) <= 0,
                    .Machine$double.eps, diag(VarCov))}
          else {
               cat("\nWARNING: Failure to solve matrix inversion of ",
                    "Approx. Hessian in LML.\n", sep="")
               cat("NOTE: Identity matrix is supplied instead.\n")
               VarCov <- diag(parm.len)}
          ### Logarithm of the Marginal Likelihood
          LML <- NA
          options(warn=-1)
          LML.test <- try(parm.len/2 * log(2*pi) + 0.5*log(det(VarCov)) +
               as.vector(Model(Modes, Data)[[1]]), silent=TRUE)
          if(!inherits(LML.test, "try-error")) LML <- LML.test[1]
          options(warn=0)
          ### Output
          LML.out <- list(LML=LML, VarCov=VarCov)
     }
     if(!is.null(theta) & !is.null(LL) & (method =="NSIS")) {
          if(!is.matrix(theta)) stop("theta must be a matrix.")
          thetacol <- ncol(theta)
          LL <- as.vector(LL)
          LLlen <- length(LL)
          if(nrow(theta) != LLlen)
               stop("The number of rows in theta differs from the ",
                    "length of LL.")
          if(LLlen < 301) {
               cat("\nWARNING: At least 301 samples are required for NSIS.\n")
               return(list(LML=NA, VarCov=NA))}
          if(thetacol > round(LLlen / 2))
               stop("At least ", round(LLlen / 2) ,
                    " stationary samples are required for",
                    thetacol," parameters for NSIS.")
          cov.prob <- 0.5
          bounds <- matrix(c(-Inf, Inf), 2, ncol(theta))
          .GetID <- function(point, center, width)
               {return(paste(floor((point - center) / width),
                    collapse=" "))}
          .PopulateHist <- function(points, heights, width)
               {
               max.point <- points[heights == max(heights), ,
                    drop=FALSE][1, ]
               center <- max.point - width / 2
               hist.env <- new.env(hash=TRUE, parent=emptyenv())
               for (i in seq(along=heights)) {
                    id <- .GetID(points[i, ], center, width)
                    if(exists(id, envir=hist.env))
                         cur.tuple <- get(id, envir=hist.env)
                    else cur.tuple <- c(Inf, 0)
                    assign(id, c(min(heights[i], cur.tuple[1]),
                         cur.tuple[2] + 1), envir=hist.env)
                    }
               return(list(env=hist.env, center=center, width=width,
                    dim=length(center)))
               }
          .Profile <- function(hist)
               {
               ids <- ls(hist$env)
               counts <- sapply(ids, function(id) get(id,
                    envir=hist$env)[2])
               names(counts) <- ids
               return(sort(counts, decreasing=TRUE))
               }
          .Lookup <- function(point, hist)
               {
               id <- .GetID(point, hist$center, hist$width)
               return(ifelse(exists(id, envir=hist$env),
                    get(id, envir=hist$env)[1], -Inf))
               }
          .SetHistNorm <- function(hist)
               {
               ids <- ls(hist$env)
               heights <- sapply(ids, function(id) get(id,
                    envir=hist$env)[1])
               max.height <- max(heights)
               hist$norm <- (log(hist$width ^ hist$dim *
                    sum(exp(heights - max.height))) + max.height)
               return(hist)
               }
          .Coverage <- function(points, hist)
               {
               heights <- apply(points, 1, function(point) .Lookup(point,
                    hist))
               return(length(heights[heights > -Inf]) / length(heights))
               }
          .Dist <- function(p1, p2)
               {
               return(max(abs(p1 - p2)))
               }
          .GetWidth <- function(theta.hist, theta.width, LL.hist,
               opt.prob=0.5)
               {
               low.point <- theta.hist[which(LL.hist == min(LL.hist))[1], ,
                    drop=FALSE]
               high.point <- theta.hist[which(LL.hist == max(LL.hist))[1], ,
                    drop=FALSE]
               minLL.dists <- apply(theta.hist, 1, function(theta) {
                    .Dist(theta, low.point)})
               maxLL.dists <- apply(theta.hist, 1, function(theta) {
                    .Dist(theta, high.point)})
               small.dist <- 0.5 * min(maxLL.dists[maxLL.dists > 0])
               big.dist <- 2 * max(minLL.dists)
               F <- function(width) {
                    hist <- .PopulateHist(theta.hist, LL.hist, width)
                    return(.Coverage(theta.width, hist) - opt.prob)
                    }
               options(warn=-1)
               solve.test <- try(uniroot(F, lower=small.dist,
                    upper=big.dist, tol=min(.05, 2/nrow(theta.width))),
                    silent=TRUE)
               options(warn=0)
               if(!inherits(solve.test, "try-error")) return(solve.test$root)
               else return(1)
               }
          .GetLimits <- function(hist)
               {
               split.ids <- strsplit(ls(hist$env), " ")
               coords <- matrix(NA, ncol=hist$dim, nrow=length(split.ids))
               for (i in seq(along=split.ids))
                    coords[i, ] <- as.integer(split.ids[[i]])
               max.coords <- apply(coords, 2, max)
               min.coords <- apply(coords, 2, min)
               result <- rbind(apply(coords, 2, min) * hist$width +
                    hist$center, (apply(coords, 2, max) + 1) * hist$width +
                    hist$center)
               rownames(result) <- c("min", "max")
               return(result)
               }
          .MargLL <- function(hist, theta.imp, LL.imp)
               {
               n <- length(LL.imp)
               samples <- rep(0, n)
               for (i in seq(length=n))
                    samples[i] <- exp(.Lookup(theta.imp[i, ],
                         hist) - LL.imp[i])
               return(list(samples=samples, mll=-log(mean(samples))))
               }
          .SplitComponents <- function(theta, LL)
               {
               tot.count <- length(LL)
               hist.count <- min(0.2 * tot.count,
                    round(sqrt(tot.count) * 2))
               width.count <- 40
               imp.count <- tot.count - hist.count - width.count
               return(list(theta.hist=theta[1:hist.count, , drop=FALSE],
                    theta.bw=theta[(hist.count + 1):(hist.count +
                         width.count), , drop=FALSE],
                    theta.imp=theta[(tot.count - imp.count + 1):tot.count, ,
                         drop=FALSE],
                    LL.hist=LL[1:hist.count],
                    LL.imp=LL[(tot.count - imp.count + 1):tot.count]))
               }
          .CheckBounds <- function(hist, bounds)
               {
               limits <- hist$limits <- .GetLimits(hist)
               if(!all((bounds[1, ] == -Inf) |
                    (bounds[1, ] < limits[1, ])))
                    stop(paste("Bounds error: histogram lower limits (",
                         paste(limits[1, ], collapse=" "),
                         ") lower than specified bounds (",
                         paste(bounds[1, ], collapse=" "), ").", sep=""))
               if(!all((bounds[2, ] == Inf) | (bounds[2, ] > limits[2, ])))
                    stop(paste("Bounds error: histogram upper limits (",
                         paste(limits[2, ], collapse=" "),
                         ") greater than specified bounds (",
                         paste(bounds[2, ], collapse=" "), ").", sep=""))
               return(hist)
               }
          options(warn=-1)
          comps <- .SplitComponents(theta, LL)
          width <- .GetWidth(comps$theta.hist, comps$theta.bw,
               comps$LL.hist, cov.prob)
          hist <- .PopulateHist(comps$theta.hist, comps$LL.hist, width)
          hist <- .SetHistNorm(hist)
          hist <- .CheckBounds(hist, bounds)
          ml.out <- .MargLL(hist, comps$theta.imp, comps$LL.imp)
          ml.out$samples <- ml.out$samples[is.finite(ml.out$samples)]
          sd.samples <- sd(ml.out$samples) / sqrt(length(ml.out$samples))
          conf.interval <- qnorm(c(.975, .025)) * sd.samples +
               mean(ml.out$samples)
          conf.interval <- ifelse(conf.interval < .Machine$double.eps,
               .Machine$double.eps, conf.interval)
          conf.interval <- -log(conf.interval)
          options(warn=0)
          LML <- ml.out$mll + hist$norm
          LML <- ifelse(is.finite(LML), LML, NA)
          ### Output
          LML.out <- list(LML=LML, VarCov=NA)
          }
     return(LML.out)
     }

#End
