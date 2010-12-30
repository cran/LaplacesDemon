###########################################################################
# Effective.Size                                                          #
#                                                                         #
# The purpose of this function is to estimate the effective sample size   #
# of a target distribution after taking autocorrelation into account.     #
# Although the code is slightly different, it is essentially the same as  #
# the effectiveSize function in the coda package.                         #
###########################################################################

Effective.Size <- function(x) 
     {
     x <- as.matrix(x)
     v0 <- order <- numeric(ncol(x))
     names(v0) <- names(order) <- colnames(x)
     z <- 1:nrow(x)
     for (i in 1:ncol(x))
          {
          lm.out <- lm(x[, i] ~ z)
          if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
               v0[i] <- 0
               order[i] <- 0
               }
          else {
               ar.out <- ar(x[, i], aic = TRUE)
               v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
               order[i] <- ar.out$order
               }
          }
     spec <- list(spec = v0, order = order)
     spec <- spec$spec
     Eff.Size <- ifelse(spec == 0, 0, nrow(x) * apply(x, 2, var)/spec)
     Eff.Size <- ifelse(Eff.Size > nrow(x), nrow(x), Eff.Size)
     return(Eff.Size)
     }
