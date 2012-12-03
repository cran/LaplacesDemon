\name{Matrices}
\alias{as.indicator.matrix}
\alias{as.inverse}
\alias{as.parm.matrix}
\alias{as.positive.definite}
\alias{as.positive.semidefinite}
\alias{as.symmetric.matrix}
\alias{Cov2Cor}
\alias{Hessian}
\alias{is.positive.definite}
\alias{is.positive.semidefinite}
\alias{is.square.matrix}
\alias{is.symmetric.matrix}
\alias{lower.triangle}
\alias{read.matrix}
\alias{tr}
\alias{upper.triangle}
\title{Matrix Utility Functions}
\description{
  These are utility functions for working with matrices.
}
\usage{
as.indicator.matrix(x)
as.inverse(x)
as.parm.matrix(x, k, parm, Data, a=-Inf, b=Inf, restrict=FALSE, chol=FALSE)
as.positive.definite(x)
as.positive.semidefinite(x)
as.symmetric.matrix(x, k=NULL)
is.positive.definite(x)
is.positive.semidefinite(x)
is.square.matrix(x)
is.symmetric.matrix(x)
Cov2Cor(Sigma)
Hessian(Model, parm, Data)
lower.triangle(x, diag=FALSE)
read.matrix(file, header=FALSE, sep="", skip=0, na.rm=FALSE)
tr(x)
upper.triangle(x, diag=FALSE)
}
\arguments{
  \item{x}{This is a matrix (though \code{as.symmetric.matrix} also
    accepts vectors).}
  \item{k}{For \code{as.parm.matrix}, this is a required argument,
    indicating the dimension of the matrix. For
    \code{as.symmetric.matrix}, this is an optional argument that
    specifies the dimension of the symmetric matrix. This applies
    only when \code{x} is a vector. It defaults to \code{NULL},
    in which case it calculates \code{k <- (-1 + sqrt(1 + 8 * length(x))) / 2}.}
  \item{diag}{Logical. If \code{TRUE}, then the elements in the main
    diagonal are also returned.}
  \item{Sigma}{This is a covariance matrix, \eqn{\Sigma}{Sigma}, and may
    be entered either as a matrix or vector.}
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{parm}{This is a vector of parameters passed to the model
    specification.}
  \item{Data}{This is the list of data passed to the model
    specification. For more information, see \code{\link{LaplacesDemon}}.}
  \item{a,b}{These optional arguments allow the elements of \code{x} to
    be bound to the interval \eqn{[a,b]}. For example, elements of a
    correlation matrix are in the interval \eqn{[-1,1]}.}
  \item{restrict}{Logical. If \code{TRUE}, then \code{x[1,1]} is
    restricted to 1. This is useful in multinomial probit, for
    example. The variable, \code{LaplacesDemonMatrix}, is
    created in a new environment, \code{LDEnv} so \code{as.parm.matrix}
    can keep track of changes from iteration to iteration.}
  \item{chol}{Logical. If \code{TRUE}, then x is an upper-triangular
    matrix.}
  \item{file}{This is the name of the file from which the data will be
    imported or read.}
  \item{header}{Logical. When \code{TRUE}, the first row of the file
    must contain names of the columns, and will be converted to the
    column names of the numeric matrix. When \code{FALSE}, the first row
    of the file contains data, not column names.}
  \item{sep}{This argument indicates a character with which it will
    separate fields when creating column vectors. For example, a
    read a comma-separated file (.csv), use \code{sep=","}.}
  \item{skip}{This is the number of rows to skip before beginning to
    read data.}
  \item{na.rm}{Logical. When \code{TRUE}, rows with missing values are
    removed from the matrix after it is read. Rather than removing
    missing values, the user may consider imputing missing values inside
    the model. Examples of imputation may be found in the accompanying
    "Examples" vignette.}
}
\details{
  The \code{as.indicator.matrix} function creates an indicator matrix
  from a vector. This function is useful for converting a discrete
  vector into a matrix in which each column represents one of the
  discrete values, and each occurence of that value in the related
  column is indicated by a one, and is otherwise filled with
  zeroes. This function is similar to the \code{class.ind} function in
  the nnet package.

  The \code{as.inverse} function returns the matrix inverse of
  \code{x}. The \code{solve} function in base R also returns the matrix
  inverse, but \code{solve} can return a matrix that is not symmetric,
  and can fail due to singularities. The \code{as.inverse} function
  tries to use the \code{solve} function to return a matrix inverse, and
  when it fails due to a singularity, \code{as.inverse} uses eigenvalue
  decomposition (in which eigenvalues below a tolerance are replaced
  with the tolerance), and coerces the result to a symmetric
  matrix. This is similar to the \code{solvcov} function in the fpc
  package.

  The \code{as.parm.matrix} function prepares a correlation, covariance,
  or precision matrix in two important ways. First,
  \code{as.parm.matrix} obtains the parameters for the matrix specified
  in the \code{x} argument by matching the name of the matrix in the
  \code{x} argument with any parameters in \code{parm}, given the
  parameter names in the \code{Data} listed in \code{parm.names}. These
  obtained parameters are organized into a matrix as the elements of the
  upper-triangular, including the diagonal. A copy is made, without the
  diagonal, and the lower-triangular is filled in, completing the
  matrix. Second, \code{as.parm.matrix} checks for
  positive-definiteness. If matrix \code{x} is positive-definite, then
  the matrix is stored as a variable called \code{LaplacesDemonMatrix}
  in a new environment called \code{LDEnv}. If matrix \code{x} is not
  positive-definite, then \code{LaplacesDemonMatrix} in \code{LDEnv} is
  sought as a replacement. If this variable exists, then it is used to
  replace the matrix. If not, then the matrix is replaced with an
  identity matrix. Back in the model specification, after using
  \code{as.parm.matrix}, it is recommended that the user also pass the
  resulting matrix back into the \code{parm} vector, so the sampler or
  algorithm knows that the elements of the matrix have changed.

  The \code{as.positive.definite} function returns the nearest
  positive-definite matrix for a matrix that is square and symmetric
  (Higham, 2002). This version is intended only for covariance and
  precision matrices, and has been optimized for speed. A more
  extensible function is \code{nearPD} in the matrixcalc package, which
  is also able to work with correlation matrices, and matrices that are
  asymmetric.

  The \code{as.positive.semidefinite} function iteratively seeks to
  return a square, symmetric matrix that is at least
  positive-semidefinite, by replacing each negative eigenvalue and
  calculating its projection. This is intended only for covariance and
  precision matrices. A similar function is \code{makePsd} in the RTAQ
  package, though it is not iterative, and returns matrices that fail a
  logical check with \code{is.positive.semidefinite}.

  The \code{as.symmetric.matrix} function accepts either a vector or
  matrix, and returns a symmetric matrix. In the case of a vector, it
  can be either all elements of the matrix, or the lower triangular. In
  the case of a \code{x} being entered as a matrix, this function
  tolerates non-finite values in one triangle (say, the lower), as long
  as the corresponding element is finite in the other (say, the upper)
  triangle.

  The \code{Cov2Cor} function converts a covariance matrix into a
  correlation matrix, and accepts the covariance matrix either in matrix
  or vector form. This function may be useful inside a model
  specification and also with converting posterior draws of the elements
  of a covariance matrix to a correlation matrix. \code{Cov2Cor} is an
  expanded form of the \code{cov2cor} function in the \code{stats}
  package, where \code{Cov2Cor} is also able to accept and return a
  vectorized matrix.

  The \code{Hessian} returns a symmetric, Hessian matrix, which is a
  matrix of second partial derivatives. The estimation of the Hessian
  matrix is approximated numerically using finite differences, as shown
  in Gelman et al. (2004, p. 313--314), though with a tolerance of
  \eqn{1e-06}. This is a slow function, and to reduce estimation time,
  the \code{Model} function is evaulated only for the upper triangular
  of this matrix, and the lower triangular is copied. This function is
  not intended to be called by the user, but is made available here.
  Automatic differentiation (AD) is an improvement over
  finite-differencing, but is currently unavailable here.

  The \code{is.positive.definite} function is a logical test of whether
  or not a matrix is positive-definite. A \eqn{k \times k}{k x k}
  symmetric matrix \eqn{\textbf{X}}{X} is positive-definite if all of
  its eigenvalues are positive (\eqn{\lambda_i > 0, i \in k}{lambda[i] >
  0, i in k}). All main-diagonal elements must be positive. The
  determinant of a positive-definite matrix is always positive, so a
  positive-definite matrix is always nonsingular. Non-symmetric,
  positive-definite matrices exist,but are not considered here.

  The \code{is.positive.semidefinite} function is a logical test of
  whether or not a matrix is positive-semidefinite. A \eqn{k x k}
  symmetric matrix \eqn{\textbf{X}}{X} is positive-semidefinite if all
  of its eigenvalues are non-negative (\eqn{\lambda_i \ge 0, i \in
  k}{lambda[i] >= 0, i ink}).

  The \code{is.square.matrix} function is a logical test of whether or
  not a matrix is square. A square matrix is a matrix with the same
  number of rows and columns, and is usually represented as a \eqn{k
  \times k}{k x k} matrix \eqn{\textbf{X}}{X}.

  The \code{is.symmetric.matrix} function is a logical test of whether
  or not a matrix is symmetric. A symmetric matrix is a square matrix
  that is equal to its transpose, \eqn{\textbf{X} = \textbf{X}^T}{X =
  t(X)}. For example, where \eqn{i} indexes rows and \eqn{j} indexes
  columns, \eqn{\textbf{X}_{i,j} = \textbf{X}_{j,i}}{X[i,j] =
  X[j,i]}. This differs from the \code{isSymmetric} function in base R
  that is inexact, using \code{all.equal}.

  The \code{lower.triangle} function returns a vector of the lower
  triangular elements of a matrix, and the diagonal is included when
  \code{diag=TRUE}.

  The \code{read.matrix} function is provided here as one of many
  convenient ways to read a numeric matrix into R. The most common
  method of storing data in R is the data frame, because it is
  versatile. For example, a data frame may contain character, factor,
  and numeric variables together. For iterative estimation, common in
  Bayesian inference, the data frame is much slower than the numeric
  matrix. For this reason, the Laplace's Demon package does not use data
  frames, and has not traditionally accepted character or factor
  data. This function is similar to the \code{read.matrix} function in
  the tseries package, with the addition of the \code{na.rm}
  argument. The \code{read.matrix} function returns a numeric matrix.

  The \code{tr} function returns the trace of a matrix. The trace of a
  matrix is the sum of the elements in the main diagonal of a square
  matrix. For example, the trace of a \eqn{k \times k}{k x k} matrix
  \eqn{\textbf{X}}{X}, is \eqn{\sum_{k=1} \textbf{X}_{k,k}}{sum(k=1)
  X[k,k]}.

  The \code{upper.triangle} function returns a vector of the lower
  triangular elements of a matrix, and the diagonal is included when
  \code{diag=TRUE}.
}
\references{
  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.
     
  Higham, N.J. (2002). "Computing the Nearest Correlation Matrix - a
  Problem from Finance". IMA Journal of Numerical Analysis, 22,
  p. 329--343.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Cov2Prec}},
  \code{\link{cov2cor}},
  \code{\link{isSymmetric}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{lower.tri}},
  \code{\link{Prec2Cov}}
  \code{\link{solve}}, and
  \code{\link{upper.tri}}.
}
\keyword{Array}
\keyword{Math}
\keyword{Matrix}
\keyword{Utility}