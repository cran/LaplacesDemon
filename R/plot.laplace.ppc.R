###########################################################################
# plot.laplace.ppc                                                        #
#                                                                         #
# The purpose of the plot.laplace.ppc function is to plot an object of    #
# class laplace.ppc.                                                      #
###########################################################################

plot.laplace.ppc <- function(x, Style=NULL, Data=NULL, Rows=NULL,
     PDF=FALSE, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(is.null(Style)) Style <- "Density"
     if(is.null(Rows)) Rows <- 1:nrow(x$yhat)
     ### Plots
     if(Style == "Covariates") {
          if(PDF == TRUE) {pdf("PPC.Plots.Covariates.pdf")
               par(mfrow=c(3,3))}
          else {par(mfrow=c(3,3), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Covariates.")
          if(is.null(Data$X) & is.null(Data$x))
               stop("X or x is required in Data.")
          if(is.null(Data$X)) co <- matrix(Data$x, length(Data$x), 1)
          else if(is.null(Data$x)) co <- Data$X
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:ncol(co)) {
               plot(co[Rows,i], temp[Rows,5], pch=16, cex=0.75,
                    ylim=c(min(temp[Rows,c(1,4:6)]),max(temp[Rows,c(1,4:6)])),
                    xlab=paste("X[,",i,"]", sep=""),
                    ylab="yhat",
                    sub="Gray dots and lines are yhat at 2.5% and 95%.")
               panel.smooth(co[Rows,i], temp[Rows,4], pch=16, cex=0.75,
                    col="gray", col.smooth="gray")
               panel.smooth(co[Rows,i], temp[Rows,6], pch=16, cex=0.75,
                    col="gray", col.smooth="gray")
               panel.smooth(co[Rows,i], temp[Rows,5], pch=16, cex=0.75)
               }
          }
     if(Style == "Covariates, Categorical DV") {
          if(PDF == TRUE) {pdf("PPC.Plots.Covariates.Cat.pdf")
               par(mfrow=c(3,3))}
          else {par(mfrow=c(3,3), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Covariates.")
          if(is.null(Data$X) & is.null(Data$x))
               stop("X or x is required in Data.")
          if(is.null(Data$X)) co <- matrix(Data$x, length(Data$x), 1)
          else if(is.null(Data$x)) co <- Data$X
          temp <- summary(x, Categorical=TRUE, Quiet=TRUE)$Summary
          ncat <- length(table(temp[,1]))
          for (i in 1:ncol(co)) {for (j in 2:(ncat+1)) {
               plot(co[Rows,i], temp[Rows,j], pch=16, cex=0.75,
                    xlab=paste("X[,",i,"]", sep=""),
                    ylab=colnames(temp)[j])
               panel.smooth(co[Rows,i], temp[Rows,j], pch=16, cex=0.75)
               }}
          }
     if(Style == "Density") {
          if(PDF == TRUE) {pdf("PPC.Plots.Density.pdf")
               par(mfrow=c(3,3))}
          else {par(mfrow=c(3,3), ask=TRUE)}
          for (j in 1:length(Rows)) {
               plot(density(x$yhat[Rows[j],]),
                    main=paste("Post. Pred. Plot of yhat[", Rows[j],
                         ",]", sep=""), xlab="Value",
                    sub="Black=Density, Red=y")
               polygon(density(x$yhat[Rows[j],]), col="black",
                    border="black")
               abline(v=x$y[Rows[j]], col="red")}
          }
     if(Style == "Fitted") {
          if(PDF == TRUE) pdf("PPC.Plots.Fitted.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          plot(temp[Rows,1], temp[Rows,5], pch=16, cex=0.75,
               ylim=c(min(temp[Rows,c(1,4:6)], na.rm=TRUE),
               max(temp[Rows,c(1,4:6)], na.rm=TRUE)),
               xlab="y", ylab="yhat", main="Fitted",
               sub="Gray dots and lines are yhat at 2.5% and 95%.")
          panel.smooth(temp[Rows,1], temp[Rows,4], pch=16, cex=0.75,
               col="gray", col.smooth="gray")
          panel.smooth(temp[Rows,1], temp[Rows,6], pch=16, cex=0.75,
               col="gray", col.smooth="gray")
          panel.smooth(temp[Rows,1], temp[Rows,5], pch=16, cex=0.75)
          }
     if(Style == "Fitted, Multivariate, C") {
          if(PDF == TRUE) {pdf("PPC.Plots.Fitted.M.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Fitted, Multivariate.")
          if(is.null(Data$Y)) stop("Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:ncol(Data$Y)) {
               plot(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[,i],
                    matrix(temp[,5], nrow(Data$Y), ncol(Data$Y))[,i],
                    pch=16, cex=0.75,
                    ylim=c(min(Data$Y[,i],
                         matrix(temp[,4], nrow(Data$Y),
                              ncol(Data$Y))[,i], na.rm=TRUE),
                         max(Data$Y[,i],
                         matrix(temp[,6], nrow(Data$Y),
                              ncol(Data$Y))[,i], na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab="yhat",
                    main="Fitted",
                    sub="Gray dots and lines are yhat at 2.5% and 95%.")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[,i],
                    matrix(temp[,4], nrow(Data$Y), ncol(Data$Y))[,i],
                    pch=16, cex=0.75, col="gray", col.smooth="gray")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[,i],
                    matrix(temp[,6], nrow(Data$Y), ncol(Data$Y))[,i],
                    pch=16, cex=0.75, col="gray", col.smooth="gray")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[,i],
                    matrix(temp[,5], nrow(Data$Y), ncol(Data$Y))[,i],
                    pch=16, cex=0.75)}
          }
     if(Style == "Fitted, Multivariate, R") {
          if(PDF == TRUE) {pdf("PPC.Plots.Fitted.M.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Fitted, Multivariate.")
          if(is.null(Data$Y)) stop("Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:nrow(Data$Y)) {
               plot(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[i,],
                    matrix(temp[,5], nrow(Data$Y), ncol(Data$Y))[i,],
                    pch=16, cex=0.75,
                    ylim=c(min(Data$Y[i,],
                         matrix(temp[,4], nrow(Data$Y),
                              ncol(Data$Y))[i,], na.rm=TRUE),
                         max(Data$Y[i,],
                         matrix(temp[,6], nrow(Data$Y),
                              ncol(Data$Y))[i,], na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab="yhat",
                    main="Fitted",
                    sub="Gray dots and lines are yhat at 2.5% and 95%.")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[i,],
                    matrix(temp[,4], nrow(Data$Y), ncol(Data$Y))[i,],
                    pch=16, cex=0.75, col="gray", col.smooth="gray")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[i,],
                    matrix(temp[,6], nrow(Data$Y), ncol(Data$Y))[i,],
                    pch=16, cex=0.75, col="gray", col.smooth="gray")
               panel.smooth(matrix(temp[,1], nrow(Data$Y), ncol(Data$Y))[i,],
                    matrix(temp[,5], nrow(Data$Y), ncol(Data$Y))[i,],
                    pch=16, cex=0.75)}
          }
     if(Style == "Predictive Quantiles") {
          if(PDF == TRUE) pdf("PPC.Plots.PQ.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          plot(temp[Rows,1], temp[Rows,7], ylim=c(0,1), pch=16, cex=0.75,
               xlab="y", ylab="PQ", main="Predictive Quantiles")
          panel.smooth(temp[Rows,1], temp[Rows,7], pch=16, cex=0.75)
          abline(h=0.025, col="gray")
          abline(h=0.975, col="gray")
          }
     if(Style == "Residuals") {
          if(PDF == TRUE) pdf("PPC.Plots.Residuals.pdf")
          par(mfrow=c(1,1))
          epsilon <- x$y - x$yhat
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          plot(epsilon.summary[2,Rows], pch=16, cex=0.75,
               ylim=c(min(epsilon.summary[,Rows], na.rm=TRUE),
                    max(epsilon.summary[,Rows], na.rm=TRUE)),
               xlab="y", ylab=expression(epsilon))
          lines(rep(0, ncol(epsilon.summary[,Rows])), col="red")
          lines(epsilon.summary[1,Rows], pch=16, cex=0.75, col="gray")
          lines(epsilon.summary[3,Rows], pch=16, cex=0.75, col="gray")
          }
     if(Style == "Residuals, Multivariate, C") {
          if(PDF == TRUE) {pdf("PPC.Plots.Residuals.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Residuals, Multivariate, C.")
          if(is.null(Data$Y))
               stop("Variable Y is required for Style=Residuals, Multivariate, C.")
          epsilon <- x$y - x$yhat
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.025 <- matrix(epsilon.summary[1,], nrow(Data$Y), ncol(Data$Y))
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data$Y), ncol(Data$Y))
          epsilon.975 <- matrix(epsilon.summary[3,], nrow(Data$Y), ncol(Data$Y))
          for (i in 1:ncol(Data$Y)) {
               plot(epsilon.500[,i], pch=16, cex=0.75,
                    ylim=c(min(epsilon.025[,i], na.rm=TRUE),
                         max(epsilon.975[,i], na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab=expression(epsilon))
               lines(rep(0, nrow(epsilon.500)), col="red")
               lines(epsilon.025[,i], pch=16, cex=0.75, col="gray")
               lines(epsilon.975[,i], pch=16, cex=0.75, col="gray")
               }
          }
     if(Style == "Residuals, Multivariate, R") {
          if(PDF == TRUE) {pdf("PPC.Plots.Residuals.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Residuals, Multivariate, C.")
          if(is.null(Data$Y))
               stop("Variable Y is required for Style=Residuals, Multivariate, C.")
          epsilon <- x$y - x$yhat
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.025 <- matrix(epsilon.summary[1,], nrow(Data$Y), ncol(Data$Y))
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data$Y), ncol(Data$Y))
          epsilon.975 <- matrix(epsilon.summary[3,], nrow(Data$Y), ncol(Data$Y))
          for (i in 1:nrow(Data$Y)) {
               plot(epsilon.500[i,], pch=16, cex=0.75,
                    ylim=c(min(epsilon.025[i,], na.rm=TRUE),
                         max(epsilon.975[i,], na.rm=TRUE)),
                    xlab=paste("Y[", i, ",]", sep=""), ylab=expression(epsilon))
               lines(rep(0, ncol(epsilon.500)), col="red")
               lines(epsilon.025[i,], pch=16, cex=0.75, col="gray")
               lines(epsilon.975[i,], pch=16, cex=0.75, col="gray")
               }
          }
     if(Style == "Space-Time by Space") {
          if(PDF == TRUE) {pdf("PPC.Plots.SpaceTime.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Space-Time by Space.")
          if(is.null(Data$longitude))
               stop("Variable longitude is required in Data.")
          if(is.null(Data$latitude))
               stop("Variable latitude is required in Data.")
          if(is.null(Data$S)) stop("Variable S is required in Data.")
          if(is.null(Data$S)) stop("Variable T is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (s in 1:Data$S) {
               plot(matrix(temp[,1], Data$S, Data$T)[s,],
                    ylim=c(min(matrix(temp[,4], Data$S, Data$T)[s,], na.rm=TRUE),
                         max(matrix(temp[,6], Data$S, Data$T)[s,]), na.rm=TRUE),
                    type="l", xlab="Time", ylab="y",
                    main=paste("Space-Time at Space s=",s," of ",
                         Data$S, sep=""),
                    sub="Actual=Black, Fit=Red, Interval=Gray")
               lines(matrix(temp[,4], Data$S, Data$T)[s,], col="gray")
               lines(matrix(temp[,6], Data$S, Data$T)[s,], col="gray")
               lines(matrix(temp[,5], Data$S, Data$T)[s,], col="red")
               }
          }
     if(Style == "Space-Time by Time") {
          if(PDF == TRUE) {pdf("PPC.Plots.SpaceTime.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Space-Time by Time.")
          if(is.null(Data$longitude))
               stop("Variable longitude is required in Data.")
          if(is.null(Data$latitude))
               stop("Variable latitude is required in Data.")
          if(is.null(Data$S)) stop("Variable S is required in Data.")
          if(is.null(Data$S)) stop("Variable T is required in Data.")
          Heat <- (1-(x$y-min(x$y, na.rm=TRUE)) /
               max(x$y-min(x$y, na.rm=TRUE), na.rm=TRUE)) * 99 + 1
          Heat <- matrix(Heat, Data$S, Data$T)
          for (t in 1:Data$T) {
               plot(Data$longitude, Data$latitude,
                    col=heat.colors(120)[Heat[,t]],
                    pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
                    main=paste("Space-Time at t=",t," of ", Data$T,
                         sep=""), sub="Red=High, Yellow=Low")}
          }
     if(Style == "Spatial") {
          if(PDF == TRUE) pdf("PPC.Plots.Spatial.pdf")
          par(mfrow=c(1,1))
          if(is.null(Data)) stop("Data is required for Style=Spatial.")
          if(is.null(Data$longitude))
               stop("Variable longitude is required in Data.")
          if(is.null(Data$latitude))
               stop("Variable latitude is required in Data.")
          heat <- (1-(x$y[Rows]-min(x$y[Rows], na.rm=TRUE)) /
               max(x$y[Rows]-min(x$y[Rows], na.rm=TRUE), na.rm=TRUE)) * 99 + 1
          plot(Data$longitude[Rows], Data$latitude[Rows],
               col=heat.colors(120)[heat],
               pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
               main="Spatial Plot", sub="Red=High, Yellow=Low")
          }
     if(Style == "Spatial Uncertainty") {
          if(PDF == TRUE) pdf("PPC.Plots.Spatial.Unc.pdf")
          par(mfrow=c(1,1))
          if(is.null(Data))
               stop("Data is required for Style=Spatial Uncertainty.")
          if(is.null(Data$longitude))
               stop("Variable longitude is required in Data.")
          if(is.null(Data$latitude))
               stop("Variable latitude is required in Data.")
          heat <- apply(x$yhat, 1, quantile, probs=c(0.025,0.975))
          heat <- heat[2,] - heat[1,]
          heat <- (1-(heat[Rows]-min(heat[Rows])) /
               max(heat[Rows]-min(heat[Rows]))) * 99 + 1
          plot(Data$longitude[Rows], Data$latitude[Rows],
               col=heat.colors(120)[heat],
               pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
               main="Spatial Uncertainty Plot",
               sub="Red=High, Yellow=Low")
          }
     if(Style == "Time-Series") {
          if(PDF == TRUE) pdf("PPC.Plots.TimeSeries.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          plot(temp[Rows,1],
               ylim=c(min(temp[Rows,c(1,4)], na.rm=TRUE),
                    max(temp[Rows,c(1,6)], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main="Plot of Fitted Time-Series",
               sub="Actual=Black, Fit=Red, Interval=Gray")
          lines(temp[Rows,4], col="gray")
          lines(temp[Rows,6], col="gray")
          lines(temp[Rows,5], col="red")
          }
     if(Style == "Time-Series, Multivariate, C") {
          if(PDF == TRUE) {pdf("PPC.Plots.TimeSeries.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Time-Series, Multivariate.")
          if(is.null(Data$Y))
               stop("Variable Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:ncol(Data$Y)) {
          plot(matrix(temp[Rows,1], nrow(Data$Y), ncol(Data$Y))[,i],
               ylim=c(min(Data$Y[,i],
                    matrix(temp[Rows,4], nrow(Data$Y),
                         ncol(Data$Y))[,i], na.rm=TRUE),
                    max(Data$Y[,i],
                    matrix(temp[Rows,6], nrow(Data$Y),
                         ncol(Data$Y))[,i], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main=paste("Time-Series ", i, " of ", ncol(Data$Y), sep=""),
               sub="Actual=Black, Fit=Red, Interval=Gray")
          lines(matrix(temp[Rows,4],nrow(Data$Y),ncol(Data$Y))[,i],
               col="gray")
          lines(matrix(temp[Rows,6],nrow(Data$Y),ncol(Data$Y))[,i],
               col="gray")
          lines(matrix(temp[Rows,5],nrow(Data$Y),ncol(Data$Y))[,i],
               col="red")}
          }
     if(Style == "Time-Series, Multivariate, R") {
          if(PDF == TRUE) {pdf("PPC.Plots.TimeSeries.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Time-Series, Multivariate.")
          if(is.null(Data$Y))
               stop("Variable Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:nrow(Data$Y)) {
          plot(matrix(temp[Rows,1], nrow(Data$Y), ncol(Data$Y))[i,],
               ylim=c(min(Data$Y[i,],
                    matrix(temp[Rows,4], nrow(Data$Y),
                         ncol(Data$Y))[i,], na.rm=TRUE),
                    max(Data$Y[i,],
                    matrix(temp[Rows,6], nrow(Data$Y),
                         ncol(Data$Y))[i,], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main=paste("Time-Series ", i, " of ", nrow(Data$Y), sep=""),
               sub="Actual=Black, Fit=Red, Interval=Gray")
          lines(matrix(temp[Rows,4],nrow(Data$Y),ncol(Data$Y))[i,],
               col="gray")
          lines(matrix(temp[Rows,6],nrow(Data$Y),ncol(Data$Y))[i,],
               col="gray")
          lines(matrix(temp[Rows,5],nrow(Data$Y),ncol(Data$Y))[i,],
               col="red")}
          }
     if(PDF == TRUE) dev.off()
     }

#End
