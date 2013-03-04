###########################################################################
# is.class                                                                #
#                                                                         #
# The purpose of the is.class functions is to provide logical tests       #
# regarding the classes of objects.                                       #
###########################################################################

is.bayesfactor <- function(x)
     {
     bayesfactor <- FALSE
     if(identical(class(x), "bayesfactor")) bayesfactor <- TRUE
     return(bayesfactor)
     }
is.bmk <- function(x)
     {
     bmk <- FALSE
     if(identical(class(x), "bmk")) bmk <- TRUE
     return(bmk)
     }
is.demonoid <- function(x)
     {
     demonoid <- FALSE
     if(identical(class(x), "demonoid")) demonoid <- TRUE
     return(demonoid)
     }
is.demonoid.hpc <- function(x)
     {
     demonoid.hpc <- FALSE
     if(identical(class(x), "demonoid.hpc")) demonoid.hpc <- TRUE
     return(demonoid.hpc)
     }
is.demonoid.ppc <- function(x)
     {
     demonoid.ppc <- FALSE
     if(identical(class(x), "demonoid.ppc")) demonoid.ppc <- TRUE
     return(demonoid.ppc)
     }
is.demonoid.val <- function(x)
     {
     demonoid.val <- FALSE
     if(identical(class(x), "demonoid.val")) demonoid.val <- TRUE
     return(demonoid.val)
     }
is.importance <- function(x)
     {
     importance <- FALSE
     if(identical(class(x), "importance")) importance <- TRUE
     return(importance)
     }
is.juxtapose <- function(x)
     {
     juxtapose <- FALSE
     if(identical(class(x), "juxtapose")) juxtapose <- TRUE
     return(juxtapose)
     }
is.laplace <- function(x)
     {
     laplace <- FALSE
     if(identical(class(x), "laplace")) laplace <- TRUE
     return(laplace)
     }
is.laplace.ppc <- function(x)
     {
     laplace.ppc <- FALSE
     if(identical(class(x), "laplace.ppc")) laplace.ppc <- TRUE
     return(laplace.ppc)
     }
is.pmc <- function(x)
     {
     pmc <- FALSE
     if(identical(class(x), "pmc")) pmc <- TRUE
     return(pmc)
     }
is.pmc.ppc <- function(x)
     {
     pmc.ppc <- FALSE
     if(identical(class(x), "pmc.ppc")) pmc.ppc <- TRUE
     return(pmc.ppc)
     }

is.pmc.val <- function(x)
     {
     pmc.val <- FALSE
     if(identical(class(x), "pmc.val")) pmc.val <- TRUE
     return(pmc.val)
     }

#End
