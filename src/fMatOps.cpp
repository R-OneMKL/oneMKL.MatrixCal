// Copyright (C) 2022             Ching-Chuan Chen
//
// This file is part of oneMKL.
//
// oneMKL.MatrixCal is free software: you can redistribute it and/or 
// modify it under the terms of the GNU General Public License as 
// published by the Free Software Foundation, either version 2 of 
// the License, or (at your option) any later version.
//
// oneMKL.MatrixCal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with oneMKL. If not, see <http://www.gnu.org/licenses/>.

#include <oneMKL.h>
#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::depends(oneMKL, RcppArmadillo)]]

//' Functions to use MKL to do the matrix calculations
//'
//' @param x,y matrices
//' @return The result matrices
//'
//' @examples
//' x <- matrix(rnorm(1e4), 100)
//' y <- matrix(rnorm(1e2), 100)
//' z <- matrix(rnorm(1e4), 100)
//' XtX <- fMatProd(t(x), x)
//' XtX2 <- fMatTransProd(x, x)
//' all.equal(XtX, XtX2) # TRUE
//' invXtX <- fMatInv(XtX)
//' fMatAdd(x, z) # x + z
//' fMatSubtract(x, z) # x - z
//' fMatSumDiffSquared(x, z) # sum((x-z)^2)
//'
//' A <- matrix(c(7,6,4,8,10,11,12,9,3,5,1,2), 3, 4)
//' A %*% fMatPseudoInv(A) # => very close to identity matrix
//' @rdname fast_matrix_ops
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatProd(const arma::mat & x, const arma::mat & y) {
  return x * y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatTransProd(const arma::mat & x, const arma::mat & y) {
  return x.t() * y;
}

//' @param is_sym_pd Whether the input matrix is symmetric/Hermitian positive definite.
//'   If the matrix is symmetric/Hermitian positive definite, enable this will be faster.
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatInv(const arma::mat & x, bool is_sym_pd = false) {
  if (is_sym_pd) {
    return arma::inv_sympd(x);
  } else {
    return arma::inv(x);
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatPseudoInv(const arma::mat & x) {
  return arma::pinv(x);
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatAdd(const arma::mat & x, const arma::mat & y) {
  return x + y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
arma::mat fMatSubtract(const arma::mat & x, const arma::mat & y) {
  return x - y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatSumDiffSquared(const arma::mat & x, const arma::mat & y) {
  return arma::sum(arma::sum(arma::square(x-y)));
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatDet(const arma::mat & x) {
  return arma::det(x);
}
