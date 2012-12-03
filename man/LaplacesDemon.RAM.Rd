\name{LaplacesDemon.RAM}
\alias{LaplacesDemon.RAM}
\title{LaplacesDemon RAM Estimate}
\description{
  This function estimates the random-access memory (RAM) required to
  update a given model and data in LaplacesDemon.

  \emph{Warning:} Unwise use of this function may crash a computer, so
  please read the details below.
}
\usage{
LaplacesDemon.RAM(Model, Data, Iterations, Thinning)
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{Data}{This is a list of Data. For more information, see
    \code{\link{LaplacesDemon}}.}
  \item{Iterations}{This is the number of iterations for which
    \code{\link{LaplacesDemon}} would update. For more information,
    see \code{\link{LaplacesDemon}}.}
  \item{Thinning}{This is the amount of thinning applied to the chains
    in \code{\link{LaplacesDemon}}.For more information, see
    \code{\link{LaplacesDemon}}.}
}
\details{
  The \code{LaplacesDemon.RAM} function uses the
  \code{\link{object.size}} function to estimate the size in MB of RAM
  required to update one chain in \code{\link{LaplacesDemon}} for a
  given model and data, and for a number of iterations and specified
  thinning. When RAM is exceeded, the computer will crash. This function
  can be useful when trying to estimate how many iterations to update a
  model without crashing the computer. However, when estimating the
  required RAM, \code{LaplacesDemon.RAM} actually creates several large
  objects, such as \code{post} (see below). If too many iterations are
  given as an argument to \code{LaplacesDemon.RAM}, for example, then it
  will crash the computer while trying to estimate the required RAM.

  The best way to use this function is as follows. First, prepare the
  model specification and list of data. Second, observe how much RAM the
  computer is using at the moment, as well as the maximum available
  RAM. The majority of the difference of these two is the amount of RAM
  the computer may dedicate to updating the model. Next, use this
  function with a small number of iterations. Note the estimated
  RAM. Increase the number of iterations, and again note the
  RAM. Continue to increase the number of iterations until, say,
  arbitrarily within 90\% of the above-mentioned difference in RAM.

  The computer operating system uses RAM, as does any other software
  running at the moment. R is currently using RAM, and other functions
  in the \code{LaplacesDemon} package, and any other package that is
  currently activated, are using RAM. There are numerous small objects
  that are not included in the returned list, that use RAM. For example,
  there may be a scalar called \code{alpha} for the acceptance
  probability, etc.

  There are potentially two larger objects that are not included, and
  the existence of the objects depends on the algorithm. First, a matrix
  used for estimating \code{\link{LML}} is used for non-adaptive MCMC
  algorithms. Second, adaptive MCMC algorithms collect a diagonal of the
  covariance matrix for each adaptation. If either is used, then the
  user should create a matrix of the appropriate dimensions and use the
  \code{\link{object.size}} function to estimate the RAM.
}
\value{
  \code{LaplacesDemon.RAM} returns a list with several components. Each
  component is an estimate in MB for an object. The list has the
  following components:
  \item{Covar}{This is the estimated size in MB of RAM required for the
    covariance matrix.}
  \item{Data}{This is the estimated size in MB of RAM required for the
    list of data.}
  \item{Deviance}{This is the estimated size in MB of RAM required for
    the deviance vector.}
  \item{Initial.Values}{This is the estimated size in MB of RAM required
    for the vector of initial values.}
  \item{Model}{This is the estimated size in MB of RAM required for the
    model specification function.}
  \item{Monitor}{This is the estimated size in MB of RAM required for
    the \eqn{N \times J}{N x J} matrix \code{Monitor}, where \eqn{N} is
    the number of thinned samples and J is the number of monitored
    variables.}
  \item{post}{This is the estimated size in MB of RAM required for the
    \eqn{N \times J}{N x J} matrix \code{post}, where \eqn{N} is the
    number of iterations and \eqn{J} is the number of initial values or
    parameters. The \code{post} matrix is not returned by
    \code{\link{LaplacesDemon}}.}
  \item{Posterior1}{This is the estimated size in MB of RAM required for
    the \eqn{N \times J}{N x J} matrix \code{Posterior1}, where \eqn{N}
    is the number of thinned samples and \eqn{J} is the number of
    initial values or parameters.}
  \item{Summary1}{This is the estimated size in MB of RAM required for
    the summary table of unthinned posterior samples of parameters,
    deviance, and monitored variables. If the chains are stationary,
    then a second table, Summary2, is also created.}
  \item{Total}{This is the estimated size in MB of RAM required in total
    to update one chain in \code{\link{LaplacesDemon}} for a given model
    and data, and for a number of iterations and specified thinning.}
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}},
  \code{\link{LML}}, and
  \code{\link{object.size}}.
}
\keyword{Memory}

