\name{indmat}
\alias{indmat}
\title{Indicator Matrix}
\description{
This function creates an indicator matrix from a vector.
}
\usage{
indmat(x)
}
\arguments{
    \item{x}{This is a vector.}
}
\value{
     \code{indmat} returns an indicator matrix in which each element of
     the matrix is zero, unless its column corresponds with the value or
     factor level position of vector \code{x}.

     This function is similar to the \code{class.ind} function in the
     \code{nnet} package.
}
\examples{

library(LaplacesDemon)
x <- c(1,2,1,3,2,4,1,5)
X <- indmat(x)
X
}
\keyword{utility}