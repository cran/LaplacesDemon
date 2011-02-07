###########################################################################
# CenterScale                                                             #
#                                                                         #
# This function centers and scales a continuous variable and provides     #
# options for binary variables. This function is very similar to Gelman's #
# rescale function in his arm package.                                    #
###########################################################################

CenterScale <- function(x, binary=NULL)
     {
     ### Initial Checks
     if(!is.numeric(x)){
          x <- as.numeric(factor(x))
          x.obs <- x[!is.na(x)]}
     x.obs <- x[!is.na(x)]
     ### Binary Variables
     if(length(unique(x.obs)) == 2){
          if(is.null(binary)){
               return((x-min(x.obs)) / (max(x.obs)-min(x.obs)))}
          else if(binary == "center0") {return (x-0.5)}
          else if(binary == "center") {return (x-mean(x.obs))}
          else if(binary == "centerscale") {
               return ((x-mean(x.obs)) / (2*sd(x.obs)))}
          }
     ### Continuous Variables
     else {return ((x-mean(x.obs)) / (2*sd(x.obs)))}
     }

#End
