## Copyright (C) 2022        Ching-Chuan Chen
##
## This file is part of oneMKL.MatrixCal.
##
## oneMKL.MatrixCal is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## oneMKL.MatrixCal is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with oneMKL.MatrixCal. If not, see <http://www.gnu.org/licenses/>.

#' `oneMKL.MatrixCal` Package
#'
#' The `oneMKL.MatrixCal` package provides an efficient and powerful tool for
#' conducting linear algebra computations in the R environment. In the package,
#' matrix operations, such as factorization, inversion, and multiplication,
#' are performed using the optimized algorithms embedded in the
#' **Intel oneAPI Math Kernel Library** (`oneMKL`). Because of this feature,
#' the computational performance can be significantly improved. To make the package
#' function work correctly with Intel `oneMKL`, users must first install the package,
#' `oneMKL`, to import Intel `oneMKL` into R. The package also supports multithreaded
#' computations, enabling users to utilize the benefit of multiple CPU cores,
#' i.e., parallel computing, to handle large-scale matrix computations.
#' The package is not only compatible with Linux, but aslo Windows operating systems,
#' which can benefit a wide range of users.
#'
#' @docType package
#' @name onemkl.MatrixCal-package
#' @useDynLib oneMKL.MatrixCal
#' @importFrom Rcpp evalCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom oneMKL onemklIncFlags onemklLibFlags
NULL
