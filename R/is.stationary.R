###########################################################################
# is.stationary                                                           #
#                                                                         #
# The purpose of the is.stationary function is to provide a logical test  #
# regarding whether or not a vector, matrix, or demonoid object is        #
# stationary. The Geweke.Diagnostic function is used.                     #
###########################################################################

is.stationary <- function(x)
  {
  if(missing(x)) stop("The x argument is required.")
  stationary <- FALSE
  if(is.vector(x)) {
       options(warn=-1)
       test <- try(as.vector(Geweke.Diagnostic(x)), silent=TRUE)
       options(warn=0)
       if(!inherits(test, "try-error") & (test > -2) & (test < 2))
            stationary <- TRUE}
  else if(is.matrix(x)) {
       options(warn=-1)
       test <- try(as.vector(Geweke.Diagnostic(x)), silent=TRUE)
       options(warn=0)
       if(!inherits(test, "try-error") & all(test > -2) & all(test < 2))
            stationary <- TRUE}
  else if(identical(class(x), "demonoid")) {
       if(x$Rec.BurnIn.Thinned < nrow(x$Posterior1)) stationary <- TRUE}
  else if(identical(class(x), "laplace")) {
       warning("x is an object of class laplace.")
       stationary <- TRUE}
  else warning("x is an unrecognized object.")
  return(stationary)
  }

#End
