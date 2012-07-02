###########################################################################
# plot.importance                                                         #
#                                                                         #
# The purpose of the plot.importance function is to plot variable         #
# importance according to the L-criterion in an object of class           #
# importance.                                                             #
###########################################################################

plot.importance <- function(x, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(!identical(class(x), "importance"))
          stop("x must be of class importance.")
     dotchart(x[,3], main="Variable Importance", xlab="L-criterion",
          pch=20)
     }

#End
