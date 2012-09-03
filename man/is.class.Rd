\name{is.class}
\alias{is.bayesfactor}
\alias{is.demonoid}
\alias{is.demonoid.hpc}
\alias{is.demonoid.ppc}
\alias{is.demonoid.val}
\alias{is.importance}
\alias{is.juxtapose}
\alias{is.laplace}
\alias{is.laplace.ppc}
\alias{is.pmc}
\alias{is.pmc.ppc}
\alias{is.pmc.val}
\title{Logical Check of Classes}
\description{
  These functions each provide a logical test of the class of an object.
}
\usage{
is.bayesfactor(x)
is.demonoid(x)
is.demonoid.hpc(x)
is.demonoid.ppc(x)
is.demonoid.val(x)
is.importance(x)
is.juxtapose(x)
is.laplace(x)
is.laplace.ppc(x)
is.pmc(x)
is.pmc.ppc(x)
is.pmc.val(x)
}
\arguments{
  \item{x}{This is an object that will be subjected to a logical test of
    its class.}
}
\details{
  Functions in Laplace's Demon often assigns a class to an output
  object. For example, after updating a model with the
  \code{\link{LaplacesDemon}} or \code{\link{LaplacesDemon.hpc}}
  function, an output object is created. The output object is of class
  \code{demonoid} or \code{demonoid.hpc}, respectively. Likewise, after
  passing a model to the \code{\link{LaplaceApproximation}} function, an
  output object is created, and it is of class \code{laplace}. The class
  of these and other objects may be logically tested.

  By assigning a class to an output object, the package is able to
  discern which other functions are appropriate for it. For example,
  after updating a model with \code{\link{LaplacesDemon}}, which creates
  an object of class \code{demonoid}, the user may desire to plot its
  output. Since it is assigned a class, the user may use the generic
  \code{plot} function, which internally selects the
  \code{\link{plot.demonoid}} function, which differs from
  \code{\link{plot.laplace}} for objects of class \code{laplace}.

  For more information on object classes, see the \code{\link{class}}
  function.
}
\value{
  The \code{is.bayesfactor} function returns a logical value indicating
  whether or not the supplied object is of class \code{bayesfactor}.
  
  The \code{is.demonoid} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid}.

  The \code{is.demonoid.hpc} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.hpc}.

  The \code{is.demonoid.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.ppc}.

  The \code{is.demonoid.val} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.val}.

  The \code{is.importance} function returns a logical value indicating
  whether or not the supplied object is of class \code{importance}.

  The \code{is.laplace} function returns a logical value indicating
  whether or not the supplied object is of class \code{laplace}.

  The \code{is.laplace.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{laplace.ppc}.

  The \code{is.pmc} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc}.

  The \code{is.pmc.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc.ppc}.

  The \code{is.pmc.val} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc.val}.
}
\references{\url{http://www.statisticat.com/laplacesdemon.html}}
\author{Byron Hall \email{statisticat@gmail.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{class}},
  \code{\link{Importance}},
  \code{\link{Juxtapose}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.laplace}},
  \code{\link{predict.pmc}}, and
  \code{\link{Validate}}.
}
\keyword{Utility}