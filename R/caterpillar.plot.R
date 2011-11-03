###########################################################################
# caterpillar.plot                                                        #
#                                                                         #
# The purpose of the caterpillar.plot function is to provide a            #
# caterpillar plot of the posterior summaries in an object of class       #
# demonoid or laplace.                                                    #
###########################################################################

caterpillar.plot <- function(x, Parms=NULL, Title=NULL)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) == "demonoid") {
          if(any(is.na(x$Summary2))) {
               x <- x$Summary1
               x.lab <- "All Samples"}
          else {
               x <- x$Summary2
               x.lab <- "Stationary Samples"}
          if(!is.null(Parms)) {
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i], rownames(x))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               if(is.numeric(Parms)) keeprows <- Parms
               temp <- x
               x <- matrix(x[keeprows,], length(keeprows), ncol(temp))
               rownames(x) <- rownames(temp)[keeprows]
               colnames(x) <- colnames(temp)
               }
          ### Setup
          x.rows <- nrow(x)
          x.lim <- c(min(x[,5]), max(x[,7]))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[,6], x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(x[i,c(5,7)], c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1 / {log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)
          }
     else if(class(x) == "laplace") {
          x.lab <- "Gaussian Samples"
          if(is.null(Parms)) {
               keeprows <- Parms <- 1:length(x$Initial.Values)}
          else {
               if(is.numeric(Parms)) keeprows <- Parms
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x$Summary))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x$Summary))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i],
                                   rownames(x$Summary))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x$Summary)))}
                         }
                    }
               }
          temp <- x$Summary
          x$Summary <- matrix(x$Summary[keeprows,], length(keeprows),
               ncol(temp))
          rownames(x$Summary) <- rownames(temp)[keeprows]
          colnames(x$Summary) <- colnames(temp)
          Modes <- x$Summary[,1]
          LB <- x$Summary[,3]
          UB <- x$Summary[,4]
          ### Setup
          x.rows <- length(Modes)
          x.lim <- c(min(LB), max(UB))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Modes
          points(Modes, x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(c(LB[i], UB[i]), c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1/{log(x.rows)/5 + 1}
          axis(2, labels=rownames(x$Summary), tick=FALSE, las=1,
               at=yy, cex.axis=cex.labels)
          }
     else {
          stop("x must be of class demonoid or laplace.")}
     }

#End
