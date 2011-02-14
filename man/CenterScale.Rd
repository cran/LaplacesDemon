\name{CenterScale}
\alias{CenterScale}
\title{Centering and Scaling}
\description{
This function centers and scales a continuous variable and provides options for binary variables.
}
\usage{
CenterScale(x, Binary=c("none","center","center0","centerscale"))
}
\arguments{
     \item{x}{This is a vector to be centered and scaled.}
     \item{Binary}{This argument indicates how binary variables will be treated, and defaults to \code{none}, which keeps the original scale, or transforms the variable to the 0-1 range, if not already there. With \code{center} in quotes, it will center the binary variable by subtracting the mean. With \code{center0} in quotes, it centers the binary variable at zero, recoding a 0 to -0.5, and a 1 to 0.5. Finally, \code{centerscale} in quotes will center and scale the binary variable, subtracting the mean and dividing by two standard deviations.}
}
\details{
Gelman (see the reference below) recommends centering and scaling continuous predictors to facilitate MCMC convergence and enable comparisons between coefficients of centered and scaled continuous predictors with coefficients of untransformed binary predictors. A continuous predictor is centered and scaled as follows: \code{x.cs <- (x - mean(x)) / (2*sd(x))}. This is an improvement over the usual practice of standardizing predictors, which is \code{x.z <- (x - mean(x)) / sd(x)}, where coefficients cannot be validly compared between binary and continuous predictors.
}
\value{
 This function returns a centered and scaled vector.
}
\references{
Gelman, A. (2008). Scaling Regression Inputs by Dividing by Two Standard
Devations. Statistics in Medicine, 27, p. 2865--2873.
}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
### See the LaplacesDemon function for an example in use.
x <- rnorm(100,10,1)
x.cs <- CenterScale(x)
}
\keyword{center, centering, predictor, scale, scaling}