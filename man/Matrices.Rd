\name{Matrices}
\alias{as.indicator.matrix}
\alias{as.inverse}
\alias{as.parm.matrix}
\alias{as.positive.definite}
\alias{as.positive.semidefinite}
\alias{as.symmetric.matrix}
\alias{is.positive.definite}
\alias{is.positive.semidefinite}
\alias{is.square.matrix}
\alias{is.symmetric.matrix}
\alias{lower.triangle}
\alias{tr}
\alias{upper.triangle}
\title{Matrix Utility Functions}
\description{
  These are utility functions for working with matrices.
}
\usage{
as.indicator.matrix(x)
as.inverse(x)
as.parm.matrix(x, k, parm, Data, a=-Inf, b=Inf, restrict=FALSE)
as.positive.definite(x)
as.positive.semidefinite(x)
as.symmetric.matrix(x, k=NULL)
is.positive.definite(x)
is.positive.semidefinite(x)
is.square.matrix(x)
is.symmetric.matrix(x)
lower.triangle(x, diag=FALSE)
tr(x)
upper.triangle(x, diag=FALSE)
}
\arguments{
  \item{x}{This is a matrix (though \code{as.symmetric.matrix} also
    accepts vectors.}
  \item{k}{For \code{as.parm.matrix}, this is a required argument,
    indicating the dimension of the matrix. For
    \code{as.symmetric.matrix}, this is an optional argument that
    specifies the dimension of the symmetric matrix. This applies
    only when \code{x} is a vector. It defaults to \code{NULL},
    in which case it calculates \code{k <- (-1 + sqrt(1 + 8 * length(x))) / 2}.}
  \item{diag}{Logical. If \code{TRUE}, then the elements in the main
    diagonal are also returned.}
  \item{parm}{This is the parameter list passed to the model
    specification.}
  \item{Data}{This is the list of data passed to the model
    specification.}
  \item{a,b}{These optional arguments allow the elements of \code{x} to
    be bound to the interval \eqn{[a,b]}. For example, elements of a
    correlation matrix are in the interval \eqn{[-1,1]}.}
  \item{restrict}{Logical. If \code{TRUE}, then \code{x[1,1]} is
    restricted to 1. This is useful in multinomial probit, for
    example. Be aware that a global variable, \code{LaplacesDemonMatrix}, is
    created so \code{as.parm.matrix} can keep track of changes from
    iteration to iteration.}
}
\details{
    The \code{as.indicator.matrix} function creates an indicator matrix
    from a vector. This function is useful for converting a discrete
    vector into a matrix in which each column represents one of the
    discrete values, and each occurence of that value in the related
    column is indicated by a one, and is otherwise filled with
    zeroes. This function is similar to the \code{class.ind} function
    in the nnet package.
    
    The \code{as.inverse} function returns the matrix inverse of
    \code{x}. The \code{solve} function in base R also returns the
    matrix inverse, but \code{solve} can return a matrix that is not
    symmetric, and can fail due to singularities. The \code{as.inverse}
    function tries to use the \code{solve} function to return a matrix
    inverse, and when it fails due to a singularity, \code{as.inverse}
    uses eigenvalue decomposition (in which eigenvalues below a
    tolerance are replaced with the tolerance), and coerces the result
    to a symmetric matrix. This is similar to the \code{solvcov}
    function in the fpc package.

    The \code{as.parm.matrix} function prepares a correlation, covariance, or
    precision matrix in two important ways. First, \code{as.parm.matrix}
    obtains the parameters for the matrix specified in the \code{x}
    argument by matching the name of the matrix in the \code{x}
    argument with any parameters in \code{parm}, given the parameter
    names in the \code{Data} listed in \code{parm.names}. These obtained
    parameters are organized into a matrix as the elements of the
    upper-triangular, including the diagonal. A copy is made, without
    the diagonal, and the lower-triangular is filled in, completing the
    matrix. Second, \code{as.parm.matrix} checks for
    positive-definiteness. If matrix \code{x} is positive-definite, then
    the matrix is stored as a global variable,
    \code{LaplacesDemonMatrix}. If matrix \code{x} is not
    positive-definite, then the global variable,
    \code{LaplacesDemonMatrix}, is sought as a replacement. If the
    global variable exists, then it is used to replace the matrix. If
    not, then the matrix is replaced with an identity matrix. Back in
    the model specification, after using \code{as.parm.matrix}, it is
    recommended that the user also pass the resulting matrix back into
    the \code{parm} vector, so the sampler or algorithm knows that
    the elements of the matrix have changed.
    
    The \code{as.positive.definite} function returns the nearest
    positive-definite matrix for a matrix that is square and symmetric,
    This is intended only for covariance and precision matrices, and has
    been optimized for speed. A more extensible function is
    \code{nearPD} in the matrixcalc package, which is also able to work
    with correlation matrices, and matrices that are asymmetric.

    The \code{as.positive.semidefinite} function iteratively seeks to
    return a square, symmetric matrix that is at least
    positive-semidefinite, by replacing each negative eigenvalue and
    calculating its projection. This is intended only for covariance
    and precision matrices. A similar function is \code{makePsd} in the
    RTAQ package, though it is not iterative, and returns matrices that fail
    a logical check with \code{is.positive.semidefinite}.

    The \code{as.symmetric.matrix} function accepts either a vector or
    matrix, and returns a symmetric matrix. In the case of a vector, it
    can be either all elements of the matrix, or the lower
    triangular. In the case of a \code{x} being entered as a matrix,
    this function tolerates non-finite values in one triangle (say, the
    lower), as long as the corresponding element is finite in the other
    (say, the upper) triangle.

    The \code{is.positive.definite} function is a logical test of
    whether or not a matrix is positive-definite. A \eqn{k \times k}{k x k}
    symmetric matrix \eqn{\textbf{X}}{X} is positive-definite if all of
    its eigenvalues are positive (\eqn{\lambda_i > 0, i \in k}{lambda[i]
    > 0, i in k}). All main-diagonal elements must be positive. The
    determinant of a positive-definite matrix is always positive, so a
    positive-definite matrix is always nonsingular. Non-symmetric,
    positive-definite matrices exist,but are not considered here.

    The \code{is.positive.semidefinite} function is a logical test of
    whether or not a matrix is positive-semidefinite. A \eqn{k x k}
    symmetric matrix \eqn{\textbf{X}}{X} is positive-semidefinite if all
    of its eigenvalues are non-negative (\eqn{\lambda_i \ge 0, i \in
    k}{lambda[i] >= 0, i in k}).

    The \code{is.square.matrix} function is a logical test of
    whether or not a matrix is square. A square matrix is a matrix with
    the same number of rows and columns, and is usually represented as a
    \eqn{k \times k}{k x k} matrix \eqn{\textbf{X}}{X}.

    The \code{is.symmetric.matrix} function is a logical test of
    whether or not a matrix is symmetric. A symmetric matrix is a square
    matrix that is equal to its transpose, \eqn{\textbf{X} =
      \textbf{X}^T}{X = t(X)}. For example, where \eqn{i} indexes rows
    and \eqn{j} indexes columns, \eqn{\textbf{X}_{i,j} =
    \textbf{X}_{j,i}}{X[i,j] = X[j,i]}. This differs from the
    \code{isSymmetric} function in base R that is inexact, using
    \code{all.equal}.

    The \code{lower.triangle} function returns a vector of the lower
    triangular elements of a matrix, and the diagonal is included when
    \code{diag=TRUE}.
        
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
  Higham, N.J. (2002). Computing the Nearest Correlation Matrix - a
    Problem from Finance. IMA Journal of Numerical Analysis, 22,
    p. 329--343.
}
\author{Byron Hall \email{laplacesdemon@statisticat.com}}
\seealso{
  \code{\link{isSymmetric}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{lower.tri}},
  \code{\link{solve}},
  \code{\link{upper.tri}}.
}
\keyword{utility}