###########################################################################
# HPD                                                                     #
#                                                                         #
# The purpose of this function is to estimate the interval of highest     #
# posterior density (HPD) of samples. Although this code is slightly      #
# different, it is essentially the same as in the coda package.           #
###########################################################################

HPD <- function(obj, prob=0.95, ...)
     {
     obj <- as.matrix(obj)
     vals <- apply(obj, 2, sort)
     if(!is.matrix(vals)) stop("obj must have nsamp > 1")
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
