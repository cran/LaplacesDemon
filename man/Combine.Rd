\name{Combine}
\alias{Combine}
\title{Combine Demonoid Objects}
\description{
This function combines objects of class \code{demonoid}, and is useful when
more than 1,000,000 iterations is required, prior to thinning, to
achieve stationarity.
}
\usage{
Combine(x, Data)
}
\arguments{
     \item{x}{This is a list of objects of class \code{demonoid}.}
     \item{Data}{This is the data, and must be identical to the data
       used to create the \code{demonoid} objects with \code{LaplacesDemon}.}
}
\details{
The \code{LaplacesDemon} algorithm limits the number of iterations to
1,000,000. Nonetheless, some complicated models may require millions of
iterations, prior to thinning, to achieve stationarity. Without using
the \code{Combine} function, Laplace's Demon (via \code{Consort}), may
suggest another million iterations to a user. After such updating,
another million may be suggested again. These suggestions may continue
to infinity, and the unnormalized logarithm of the joint posterior
density, \code{LP}, may seem to oscillate, sometimes improving and
getting higher, and getting lower during other updates.

In these cases, which usually involve complicated models with high
autocorrelation in the chains, the user should save (meaning, not
overwrite) each object of class \code{demonoid}, place multiple objects
into a list, and use the \code{Combine} function to combine these
objects.

For example, suppose a user names the object Fit, as in the
\code{LaplacesDemon} example. Now, rather than overwriting object Fit,
object Fit is renamed, after updating a million iterations, to Fit1. As
suggested by \code{Consort}, another million iterations are used, but
now to create object Fit2. Further suppose this user specified
\code{Thinning=1000} in \code{LaplacesDemon}, meaning that the million
iterations are thinned by 1,000, so only 1,000 iterations are retained
in each object, Fit1 and Fit2. In this case, \code{Combine} combines the
information in Fit1 and Fit2, and returns an object the user names
Fit3. Fit3 has only 1,000 iterations, which is the result of appending
the iterations in Fit1 and Fit2, and thinning by 2. If 2,000,000
iterations were updated from the beginning (which cannot be, due to
the limit of 1,000,000 iterations), and were thinned by 2,000, then the
same information exists now in Fit3. The \code{Consort} function can now
be applied to Fit3, to see if stationarity is found. If not, more
objects of class \code{demonoid} can be collected and combined.

This is faster than updating millions of iterations in one function call.
}
\value{
 This function returns an object of class \code{demonoid}. For more
 information on an object of class \code{demonoid}, see the
 \code{LaplacesDemon} function.
}
\seealso{\code{\link{LaplacesDemon}}}

