###########################################################################
# CenterScale                                                             #
#                                                                         #
# This function centers and scales a continuous variable and provides     #
# options for binary variables. This function is very similar to Gelman's #
# rescale function in his arm package.                                    #
###########################################################################

CenterScale <- function(x, Binary=c("none", "center", "center0",
     "centerscale"))
     {
     ### Initial Checks
     if(!is.numeric(x)){
          x <- as.numeric(factor(x))
          x.obs <- x[!is.na(x)]}
     x.obs <- x[!is.na(x)]
     if(is.null(Binary)) Binary="none"
     ### Binary Variables
     if(length(unique(x.obs)) == 2){
          if(Binary == "none"){
               return((x-min(x.obs)) / (max(x.obs)-min(x.obs)))}
          else if(Binary == "center") {return (x-mean(x.obs))}
          else if(Binary == "center0") {
               x <- (x-min(x.obs)) / (max(x.obs)-min(x.obs))
               return (x-0.5)}
          else if(Binary == "centerscale") {
               return ((x-mean(x.obs)) / (2*sd(x.obs)))}
          }
     ### Continuous Variables
     else {return ((x-mean(x.obs)) / (2*sd(x.obs)))}
     }

#End
