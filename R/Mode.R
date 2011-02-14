###########################################################################
# Mode                                                                    #
#                                                                         #
# The purpose of this function is to return the mode of a vector.         #
###########################################################################

Mode <- function(x)
     {
     x <- as.vector(x)
     ### Discrete
     if((length(unique(x)) < length(x)) & (length(unique(x)) <= 20)) {
     Mode <- as.vector(which.max(table(x)))}
     ### Continuous (using kernel density)
     else {
          x <- as.vector(as.numeric(as.character(x)))
          kde <- density(x)
          Mode <- kde$x[kde$y==max(kde$y)]}
     return(Mode)
     }

#End
