% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fast_matrix_extra_ops}
\alias{fast_matrix_extra_ops}
\alias{fMatAdd}
\alias{fMatSubtract}
\alias{fMatElementWiseProduct}
\alias{fMatElementWiseDivide}
\alias{fMatRowSum}
\alias{fMatColSum}
\alias{fMatRowMin}
\alias{fMatColMin}
\alias{fMatRowMax}
\alias{fMatColMax}
\alias{fMatSumDiffSquared}
\title{Functions that use oneMKL for fast matrix calculations through RcppEigen}
\usage{
fMatAdd(X, Y)

fMatSubtract(X, Y)

fMatElementWiseProduct(X, Y)

fMatElementWiseDivide(X, Y)

fMatRowSum(X)

fMatColSum(X)

fMatRowMin(X)

fMatColMin(X)

fMatRowMax(X)

fMatColMax(X)

fMatSumDiffSquared(X, Y)
}
\arguments{
\item{X, Y}{The input matrices 'X' and 'Y'.}
}
\value{
The corresponding results.
}
\description{
\describe{
\item{\strong{fMatAdd}}{This function returns the sum of matrices `X` and `Y`, i.e., `X + Y`}
\item{\strong{fMatSubtract}}{This function returns the result of the matrix `X` minus the matrix `Y`, namely, `X - Y`.}
\item{\strong{fMatElementWiseProduct}}{This function returns the result of the matrix `X` element-wisely product the matrix `Y`, i.e. [Xij * Yij] for all i, j.}
\item{\strong{fMatElementWiseDivide}}{This function returns the result of the matrix `X` element-wisely divided by the matrix `Y`, i.e. [Xij / Yij] for all i, j.}
\item{\strong{fMatRowSum}}{This function returns the sum of each row.}
\item{\strong{fMatColSum}}{This function returns the sum of each column.}
\item{\strong{fMatRowMin}}{This function returns the minimum of each row.}
\item{\strong{fMatRowMax}}{This function returns the maximum of each row.}
\item{\strong{fMatColMin}}{This function returns the minimum of each column.}
\item{\strong{fMatColMax}}{This function returns the maximum of each column.}
\item{\strong{fMatSumDiffSquared}}{The function returns the result of square of `X` minus `Y`, i.e., `(X-Y)^2`, where `X` and `Y` are matrices.}
}
}
\examples{
x <- matrix(rnorm(1e4), 100)
y <- matrix(rnorm(1e2), 100)
z <- matrix(rnorm(1e4), 100)

fMatAdd(x, z) # x + z
fMatSubtract(x, z) # x - z

fMatElementWiseProduct(x, z) # x * z
fMatElementWiseDivide(x, z) # x / z

fMatColSum(x) # colSums(x)
fMatRowSum(x) # rowSums(x)
fMatRowMin(x) # apply(x, 1, min)
fMatRowMax(x) # apply(x, 1, max)
fMatColMin(x) # apply(x, 2, min)
fMatColMax(x) # apply(x, 2, max)

fMatSumDiffSquared(x, z) # sum((x-z)^2)
}
