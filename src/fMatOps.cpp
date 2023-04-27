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

// [[Rcpp::depends(oneMKL)]]

//' Functions that use oneMKL for fast matrix calculations
//'
//' @param X,Y matrices
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
//' fMatRowSum(x) # rowSums(x)
//' fMatRowMin(x) # apply(x, 1, min)
//' fMatRowMax(x) # apply(x, 1, max)
//' fMatColMin(x) # apply(x, 2, min)
//' fMatColMax(x) # apply(x, 2, max)
//' @rdname fast_matrix_ops
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatProd(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y
) {
  return X * Y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatTransProd(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y
) {
  return X.transpose() * Y;
}

//' @param is_invertible specify whether to enable faster computation of the linear model solution
//'   by disabling the use of rcond, iterative refinement, and equilibration.
//' @param is_sym_pd specific whether the input matrix is symmetric/Hermitian positive definite.
//'   Enabling this option can result in faster computation if the matrix satisfies these properties.
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatSolve(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y,
    bool is_sym_pd = false,
    bool is_invertible = false
) {
  if (is_sym_pd) {
    return X.llt().solve(Y);
  } else if (is_invertible) {
    return X.partialPivLu().solve(Y);
  } else {
    return X.fullPivLu().solve(Y);
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatInv(const Eigen::Map<Eigen::MatrixXd> X, bool is_sym_pd = false) {
  if (is_sym_pd) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X.rows(), X.cols());
    return X.llt().solve(I);
  } else {
    return X.inverse();
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatAdd(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y
) {
  return X + Y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatSubtract(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y
) {
  return X - Y;
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatRowSum(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.rowwise().sum();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatColSum(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.colwise().sum();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatRowMin(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.rowwise().minCoeff();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatColMin(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.colwise().minCoeff();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatRowMax(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.rowwise().maxCoeff();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatColMax(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.colwise().maxCoeff();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatDet(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.determinant();
}
