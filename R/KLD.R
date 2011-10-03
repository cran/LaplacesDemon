###########################################################################
# Kullback-Leibler Divergence (KLD)                                       #
#                                                                         #
# The purpose of the KLD() function is to calculate the Kullback-Leibler  #
# divergences between two probability distributions, p(x) and p(y).       #
###########################################################################

KLD <- function (px, py, base=exp(1))
     {
     ### Initial Checks
     if(is.matrix(px) && ncol(px) == 2) px <- px[, 2]
     if(is.matrix(py) && ncol(py) == 2) py <- py[, 2]
     n1 <- length(px)
     n2 <- length(py)
     if(n1 != n2) stop("px and py must have the same length.")
     if(any(is.na(px)) | any(is.na(py))) {
          KLD.px.py <- KLD.py.px <- NA
          warning("Missing values found; KLD is set to missing.",
               call.=FALSE)
          }
     else {
          if(any(px < 0, na.rm=TRUE) | any(py < 0, na.rm=TRUE))
               stop("Negative density found in px or py.")
          if(sum(px) == 0) {
               warning("WARNING: sum(px) = 0.", call.=FALSE)
               px[px == 0] <- 1e-07
               }
          if(sum(py) == 0) {
               warning("WARNING: sum(py) = 0.", call.=FALSE)
               py[py == 0] <- 1e-07}
          }
     ### Normalize
     px <- px / sum(px)
     py <- py / sum(py)
     ### Kullback-Leibler Calculations
     KLD.px.py <- px * log(px / py, base=base)
     KLD.py.px <- py * log(py / px, base=base)
     sum.KLD.px.py <- sum(KLD.px.py)
     sum.KLD.py.px <- sum(KLD.py.px)
     mean.KLD <- (KLD.px.py + KLD.py.px) / 2
     mean.sum.KLD <- (sum.KLD.px.py + sum.KLD.py.px) / 2
     ### Output
     out <- list(KLD.px.py=KLD.px.py, #KLD[i](p(x[i]) || p(y[i]))
          KLD.py.px=KLD.py.px, #KLD[i](p(y[i]) || p(x[i]))
          mean.KLD=mean.KLD,
          sum.KLD.px.py=sum.KLD.px.py, #KLD(p(x) || p(y))
          sum.KLD.py.px=sum.KLD.py.px, #KLD(p(y) || p(x))
          mean.sum.KLD=mean.sum.KLD,
          intrinsic.discrepancy=min(sum.KLD.px.py, sum.KLD.py.px))
     return(out)
     }

#End
