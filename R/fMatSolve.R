## Copyright (C) 2022-2023        Ching-Chuan Chen
##
## This file is part of oneMKL.MatrixCal.
##
## oneMKL.MatrixCal is free software: you can redistribute it and/or 
## modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation, either version 2 of 
## the License, or (at your option) any later version.
##
## oneMKL.MatrixCal is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with oneMKLUtil. If not, see <http://www.gnu.org/licenses/>.

#' Replaced base::solve Function
#'
#' @param a,b,tol,LINPACK Please refer to \link[base]{solve}.
#'
#' @export
fMatSolve <- function(a, b, tol = .Machine$double.eps, LINPACK = FALSE, ...) {
  if(is.complex(a) || (!missing(b) && is.complex(b))) {
    a <- as.matrix(a)
    if (missing(b)) {
      b <- diag(1.0+0.0i, nrow(a))
    }
    return(mkl_cmpl_solve(a, as.matrix(b)))
  }

  if (is.qr(a)) {
    warning("solve.default called with a \"qr\" object: use 'qr.solve'")
    return(solve.qr(a, b, tol))
  }

  a <- as.matrix(a)
  if (missing(b)) {
    b <- diag(1.0, nrow(a))
  }

  return(mkl_real_solve(a, as.matrix(b), tol))
}
