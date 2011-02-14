\name{Mode}
\alias{Mode}
\title{The Mode of a vector}
\description{
The mode is a measure of central tendency. It is the value that occurs
most frequently, or in a continuous probability distribution, it is the
value with the most density.
}
\usage{
Mode(x)
}
\arguments{
     \item{x}{This is a vector in which a single mode will be
       sought.}
}
\details{
When vector \code{x} has a number of unique values that is less than its
length and less than or equal to 20, the \code{Mode} function treats
vector \code{x} as discrete and returns the most frequent value. If
multiple modes exist, this function will return only the first one
found.

When vector \code{x} has a number of unique values that is greater than
20 or equal to its length, the \code{Mode} function treats vector{x} as
a continuous variable, and returns the value of \code{x} associated with
the most kernel density. If multiple modes exist, this function will
return only the first one found.
}
\value{
 \code{Mode} returns the value that occurs most frequently, or in a
 continuous probability distribution, the value with the most kernel
 density.
}
\seealso{\code{\link{LaplacesDemon}}}
\examples{
x <- rnorm(100,0,1)
Mode(x)
}
\keyword{mode}