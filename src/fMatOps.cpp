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
#include "utils.h"

// [[Rcpp::depends(oneMKL)]]

//' Functions that use oneMKL for fast matrix calculations through RcppEigen
//'
//' \describe{
//' \item{\strong{fMatProd}}{This function returns the multiplication of matrices `X` and `Y`, i.e., `XY`.}
//' \item{\strong{fMatTransProd}}{This function returns the product of the transpose of the matrix `X`
//'  and the matrix `Y`, i.e., `X^T Y`.}
//' \item{\strong{fMatSolve}}{This function returns the solution of a linear system `AX=Y`.
//'  If the matrix `X` is symmetric positive definite, Cholesky decomposition
//'  will be used for better computational performance.
//'  If the matrix `X` is invertible, the LU decomposition
//'  will be used for better computational performance.}
//' \item{\strong{fMatInv}}{This function returns the inverse of the matrix `X`, i.e., `X^(-1)`.
//'  If the matrix `X` is symmetric positive definite, Cholesky decomposition
//'  will be used for better computational performance.}
//' \item{\strong{fMatPseudoInv}}{This function returns the pseudo-inverse (also called generalized inverse or g-inverse) of the matrix 'X'.}
//' \item{\strong{fMatLeastSquare}}{This function returns the solution of least square through QR decomposition.
//' If `stable` = TRUE, a more stable but slower algorithm will be applied.}
//' \item{\strong{fMatDet}}{This function returns the determinant of the matrix `X`.}
//' \item{\strong{fMatRank}}{This function returns the rank of the matrix `X`.}
//' \item{\strong{fMatRCond}}{This function returns the reciprocal condition number of the matrix `X`.}
//' }
//'
//' @param X,Y The input matrices 'X' and 'Y'.
//' @param is_X_symmetric A logical variable indicating whether the input matrix `X` is symmetric.
//'  Better computational performance is expected if the matrix is symmetric.
//' @return The corresponding results.
//'
//' @examples
//' x <- matrix(rnorm(1e4), 100)
//' y <- matrix(rnorm(1e2), 100)
//' z <- matrix(rnorm(1e4), 100)
//' XtX <- fMatProd(t(x), x)
//' XtX2 <- fMatTransProd(x, x)
//' all.equal(XtX, XtX2) # TRUE
//'
//' invXtX <- fMatInv(XtX)
//' fMatSolve(XtX, fMatTransProd(x, y)) # linear regression coefficients
//'
//' fMatDet(x)
//' fMatRank(x)
//' fMatRCond(x)
//' @rdname fast_matrix_ops
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatProd(SEXP X, SEXP Y, bool is_X_symmetric = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if (Rf_ncols(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of columns of X");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (Rf_ncols(Y) == 1) {
    // to have the best performance
    Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * b;
    } else {
      return XMtd * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return XMtd * Z;
    }
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatTransProd(SEXP X, SEXP Y, bool is_X_symmetric = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of rows of X");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  int *ydims;
  ydims = INTEGER(Rf_coerceVector(Rf_getAttrib(Y, R_DimSymbol), INTSXP));
  if (ydims[1] == 1) {
    // to have the best performance
    Eigen::Map<Eigen::VectorXd> b = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * b;
    } else {
      return XMtd.transpose() * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    if (is_X_symmetric) {
      return XMtd.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return XMtd.transpose() * Z;
    }
  }
}

//' @param is_invertible A logical variable indicating whether the input matrix `X` is invertible.
//'  Better computational performance is expected if the matrix is invertible for `fMatSolve`.
//' @param is_sym_pd A logical variable indicating whether the input matrix `X` is symmetric positive definitive.
//'  Better computational performance is expected if the matrix is symmetric positive definitive for `fMatSolve` and `fMatInv`.
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatSolve(
    SEXP X, SEXP Y,
    bool is_sym_pd = false,
    bool is_invertible = false
) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  if (Rf_nrows(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of columns/rows of X");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
  if (is_sym_pd) {
    return XMtd.llt().solve(YMtd);
  } else if (is_invertible) {
    return XMtd.partialPivLu().solve(YMtd);
  } else {
    return XMtd.householderQr().solve(YMtd);
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatInv(SEXP X, bool is_sym_pd = false) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (is_sym_pd) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(XMtd.rows(), XMtd.cols());
    return XMtd.llt().solve(I);
  } else {
    return XMtd.inverse();
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatPseudoInv(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.completeOrthogonalDecomposition().pseudoInverse();
}

//' @param stable A logical variable indicating whether to use a more stable
//'  but slower algorithm for `fMatLeastSquare`.
//' @param is_X_full_rank A logical variable indicating whether the input matrix 'X' is full-rank.
//'  If false, the pseudo inverse will be applied for `fMatLeastSquare`.
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatLeastSquare(
    SEXP X, SEXP Y,
    bool stable = false,
    bool is_X_full_rank = true
) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_nrows(Y)) {
    Rcpp::stop("The number of rows of Y must be equal to the number of rows of X");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
  if (is_X_full_rank) {
    if (stable) {
      return XMtd.householderQr().solve(YMtd);
    } else {
      return XMtd.colPivHouseholderQr().solve(YMtd);
    }
  } else {
    return XMtd.completeOrthogonalDecomposition().solve(YMtd);
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatDet(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.determinant();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatRank(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.colPivHouseholderQr().rank();
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatRCond(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.partialPivLu().rcond();
}
