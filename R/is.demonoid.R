###########################################################################
# is.demonoid; is.laplace                                                 #
#                                                                         #
# The purpose of the is.demonoid and is.laplace functions is to provide   #
# logical tests regarding the classes of objects.                         #
###########################################################################

is.demonoid <- function(x)
     {
     demonoid <- FALSE
     if(class(x) == "demonoid") demonoid <- TRUE
     return(demonoid)
     }
is.laplace <- function(x)
     {
     laplace <- FALSE
     if(class(x) == "laplace") laplace <- TRUE
     return(laplace)
     }

#End
