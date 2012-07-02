###########################################################################
# is.proper                                                               #
#                                                                         #
# The purpose of the is.proper function is to provide a logical check of  #
# whether or not a probability distribution is proper, meaning whether or #
# not it integrates to one.                                               #
###########################################################################

is.proper <- function(f, a, b, tol=1e-5)
       {
       ### Initial Checks
       if(!is.function(f) & (class(f) != "demonoid") &
            (class(f) != "laplace"))
            stop("f is not a function or object of class demonoid or laplace.")
       ### Propriety
       propriety <- FALSE
       if(is.function(f)) {
            if(a >= b) stop("a >= b.")
            area <- integrate(f,a,b)$value
            if((area >= (1-tol)) & (area <= (1+tol))) propriety <- TRUE}
       else if(is.finite(f$LML)) propriety <- TRUE
       return(propriety)
       }

#End
