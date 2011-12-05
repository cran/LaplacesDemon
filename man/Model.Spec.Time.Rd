\name{Model.Specification.Time}
\alias{Model.Spec.Time}
\title{Model Specification Time}
\description{
  The \code{Model.Spec.Time} function returns the time in minutes that it
  took to evaluate a model specification a number of times, as well as
  the evaluations per minute, and componentwise iterations per minute.
}
\usage{
Model.Spec.Time(Model, Initial.Values, Data, n=1000)
}
\arguments{
  \item{Model}{This requried argument is a model specification
    function. For more information, see \code{\link{LaplacesDemon}}.}
  \item{Initial.Values}{This required argument is a vector of initial
    values for the parameters.}
  \item{Data}{This required argument is a list of data. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{n}{This is the number of evaluations of the model specification,
    and accuracy increases with \code{n}.}
}
\details{
  The largest single factor to affect the run-time of an algorithm,
  whether it is with \code{\link{LaplaceApproximation}} or
  \code{\link{LaplacesDemon}}, is the time that it takes to evaluate the
  model specifiation. This has also been observed in Rosenthal (2007).
  It is highly recommended that a user of the \code{LaplacesDemon}
  package should attempt to reduce the run-time of the model
  specification, usually by testing alternate forms of code for speed.

  A model specification function that is optimized for speed and
  involves many records may result in a model update run-time that is
  considerably less than other popular MCMC-based software algorithms
  that loop through records, even when those algorithms are coded in
  \code{C} or other fast languages. For a comparison, see the
  ``Laplace's Demon Tutorial'' vignette.

  However, if a model specification function in the \code{LaplacesDemon}
  package is not fully vectorized (contains \code{for} loops and
  \code{apply} functions), then run-time will typically be slower than
  other popular MCMC-based software algorithms.

  Examples in the \code{LaplacesDemon} package are not necessarily
  optimized for speed, which is often compromised for clarity of
  programming. For example, most examples in the ``Examples'' vignette
  use the \code{grep} function to find the position of a parameter based
  on an associated parameter name. This makes the model specification
  function much more readable, and also makes it easy to modify by
  changing the number and position of parameters, such as when a user
  experiments with the model form, or changes it from the example data
  to actual data. However, the \code{grep} function is much slower than
  specifying the exact location. There are many other instances in the
  examples that demonstrate user-friendly code and are not optimized for
  speed.

  When speed is of paramount importance, a user is encouraged to
  experiment with writing the model specification function in another
  language such as \code{C} or \code{Fortran}, and calling the new
  function from within the \code{Model} function with the \code{.C} or
  \code{.Fortran} function. This experimentation has not yet been done,
  but should generally be worth the effort only when the model
  specification is vectorized in the faster language as well. One
  suggestion may be to browse for R packages on CRAN that enable a user
  to speed up R code by incorporating C and other languages.

  Another use for \code{Model.Spec.Time} is to allow the user to make an
  informed decision about which MCMC algorithm to select, given the
  speed of their model specification. For example, the Adaptive
  Metropolis-within-Gibbs (AMWG) of Roberts and Rosenthal (2009) is
  currently the adaptive algorithm of choice in general, but this choice
  is conditional on run-time. While the other MCMC algorithms in
  \code{LaplacesDemon} evaluate the model specification function once
  per iteration, componentwise algorithms (AMWG and MWG) evaluate it
  once per parameter per iteration, significantly increasing run-time
  per iteration in large models. The \code{Model.Spec.Time} may forewarn
  the user of the associated run-time, and if it should be decided to go
  with an alternate MCMC algorithm, such as AMM, in which each element
  of its covariance matrix must stabilize for the chains to become
  stationary. AMM, for example, will require many more iterations of
  burn-in (for more information, see the \code{\link{burnin}} function),
  but with numerous iterations, allows more thinning. A general
  recommendation may be to use AMWG when
  \code{Componentwise.Iters.per.Minute} >= 1000, but this is subjective
  and best determined by each user for each model. This does not include
  the time of the MCMC algorithm, which includes pre- and post-processing
  time (checking stationarity, summarizing, estimating LML, etc.), and
  therefore changes with the number of iterations.

  Following are a few common suggestions for increasing the speed of R
  code in the model specification function. There are often exceptions
  to these suggestions, so case-by-case experimentation is also suggested.

  \itemize{
    \item Replace exponents with multiplied terms, such as \code{x^2}
    with \code{x*x}.
    \item Replace \code{mean(x)} with \code{sum(x)/length(x)}.
    \item Replace parentheses (when possible) with curly brackets.
    \item Specify parameters directly, such as \code{parm[53]}, rather
    than with \code{grep}.
    \item Vectorize functions and eliminate \code{apply} and \code{for}
    functions. There are often specialized functions. For example,
    \code{max.col(X)} is faster than \code{apply(X, 1, which.max)}.
  }
}
\value{
  The \code{Model.Spec.Time} function returns a list with three
  components:

  \item{Time}{This is the time in minutes to evaluate the model
    specification \code{n} times.}
  \item{Evals.per.Minute}{This is the number of evaluations of the
    model specification per minute: \code{n} / \code{Time}. This is also
    the number of iterations per minute in an algorithm that is not
    componentwise, where one evaluation occurs per iteration.}
  \item{Componentwise.Iters.per.Minute}{This is the number of iterations
    per minute in a componentwise MCMC algorithm, such as AMWG or
    MWG. It is the evaluations per minute divided by the number of
    parameters, since an evaluation must occur for each parameter, for
    each iteration.}
}
\references{
  Hall, B. (2011). "Laplace's Demon", STATISTICAT, LLC.
  URL=\url{http://www.statisticat.com/laplacesdemon.html}

  Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive
  MCMC". Computational Statistics and Data Analysis, 18, p. 349--367.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{.C}},
  \code{\link{.Fortran}},
  \code{\link{apply}},
  \code{\link{grep}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{max.col}},
  \code{\link{system.time}}.
}
\keyword{utility}
