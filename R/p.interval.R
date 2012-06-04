###########################################################################
# p.interval                                                              #
#                                                                         #
# The purpose of the p.interval function is to estimate either a          #
# quantile-based probability interval, uni-modal highest posterior        #
# density (HPD) interval, or multi-modal HPD intervals of posterior       #
# samples. The code for the uni-modal HPD interval is similar to code in  #
# the coda package, but this is designed to work with matrices or objects #
# of class demonoid or laplace.                                           #
###########################################################################

p.interval <- function(obj, HPD=TRUE, MM=TRUE, prob=0.95, ...)
     {
     ### Initial Checks
     if(missing(obj)) stop("The obj argument is required.")
     if(any(!is.finite(obj)))
          stop("The obj argument must contain finite values.")
     if(length(unique(as.vector(obj))) <= 1)
          stop("The obj argument has insufficient length.")
     if(length(prob) > 1)
          stop("The prob argument must be a scalar.")
     if((prob <= 0) | (prob >= 1))
          stop("The prob argument must be in the interval (0,1).")
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
     else if(class(obj) == "laplace") {
          if(any(is.na(obj$Posterior)))
               stop("Posterior samples do not exist in obj.")
          obj <- obj$Posterior
          }
     else obj <- as.matrix(obj)
     ### Probability Intervals
     if(HPD == TRUE) {
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
          attr(ans, "Probability.Interval") <- gap / nsamp
          }
     if(HPD == FALSE) {
          ans <- apply(obj, 2, quantile,
               probs=c((1-prob)/2,1-((1-prob)/2)))
          ans <- as.matrix(t(ans))
          colnames(ans) <- c("Lower","Upper")
          attr(ans, "Probability.Interval") <- prob}
     if(MM == TRUE) {
          mm <- apply(obj, 2, is.multimodal)
          if(any(mm)) {
               cat("\n\nPotentially multi-modal column vectors:\n",
                    which(MM),"\n")
               for (m in which(mm)) {
                    kde <- density(vals[,m])
                    dens <- approx(kde$x, kde$y, vals[,m])$y
                    dens.ind <- dens >= as.vector(quantile(dens,
                         probs=1-prob)) * 1
                    ints <- ""
                    for (i in 1:nrow(vals)) {
                         if((i == 1) & (dens.ind[i] == 1)) 
                              ints <- paste("(",round(vals[i,m],3),",",sep="")
                         if(i > 1) {
                              if((dens.ind[i] == 0) & (dens.ind[i-1] == 1)) 
                                   ints <- paste(ints,round(vals[i-1,m],3),")",sep="")
                              if((dens.ind[i] == 1) & (dens.ind[i-1] == 0)) 
                                   ints <- paste(ints," (",round(vals[i,m],3),",",sep="")}
                         }
                    if((dens.ind[i] == 1) & (dens.ind[i-1] == 1)) {
                         ints <- paste(ints,round(vals[i,m],3),")",sep="")}
                    cat("\nColumn", m, "multimodal intervals:", ints, "\n")}
               }
          }
     return(ans)
     }

#End
