###########################################################################
# HPD                                                                     #
#                                                                         #
# The purpose of the HPD function is to estimate the interval of highest  #
# posterior density (HPD) of samples. This function is similar to code in #
# the coda package, but this is designed to work with objects of class    #
# demonoid.                                                               #
###########################################################################

HPD <- function(obj, prob=0.95, ...)
     {
     if(class(obj) == "demonoid") {
          thin <- obj$Rec.BurnIn.Thinned
          n <- nrow(obj$Posterior1)
          if(thin < nrow(obj$Posterior1))
               x <- cbind(obj$Posterior1[thin:n,], obj$Monitor[thin:n,])
          if(thin >= nrow(obj$Posterior1))
               x <- cbind(obj$Posterior1, obj$Monitor)
          colnames(x) <- c(colnames(obj$Posterior1),
               colnames(obj$Monitor))
          obj <- x
          }
     else obj <- as.matrix(obj)
     vals <- apply(obj, 2, sort)
     if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
     nsamp <- nrow(vals)
     npar <- ncol(vals)
     gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
     init <- 1:(nsamp - gap)
     inds <- apply(vals[init + gap, , drop=FALSE] -
          vals[init, , drop=FALSE], 2, which.min)
     ans <- cbind(vals[cbind(inds, 1:npar)],
          vals[cbind(inds + gap, 1:npar)])
     dimnames(ans) <- list(colnames(obj), c("Lower", "Upper"))
     attr(ans, "Probability") <- gap / nsamp
     return(ans)
     }

#End
