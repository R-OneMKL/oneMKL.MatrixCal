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
//' \item{\strong{fMatAdd}}{This function returns the sum of matrices `X` and `Y`, i.e., `X + Y`}
//' \item{\strong{fMatSubtract}}{This function returns the result of the matrix `X` minus the matrix `Y`, namely, `X - Y`.}
//' \item{\strong{fMatElementWiseProduct}}{This function returns the result of the matrix `X` element-wisely product the matrix `Y`, i.e. [Xij * Yij] for all i, j.}
//' \item{\strong{fMatElementWiseDivide}}{This function returns the result of the matrix `X` element-wisely divided by the matrix `Y`, i.e. [Xij / Yij] for all i, j.}
//' \item{\strong{fMatRowSum}}{This function returns the sum of each row.}
//' \item{\strong{fMatColSum}}{This function returns the sum of each column.}
//' \item{\strong{fMatRowMin}}{This function returns the minimum of each row.}
//' \item{\strong{fMatRowMax}}{This function returns the maximum of each row.}
//' \item{\strong{fMatColMin}}{This function returns the minimum of each column.}
//' \item{\strong{fMatColMax}}{This function returns the maximum of each column.}
//' \item{\strong{fMatSumDiffSquared}}{The function returns the result of square of `X` minus `Y`, i.e., `(X-Y)^2`, where `X` and `Y` are matrices.}
//' }
//'
//' @param X,Y The input matrices 'X' and 'Y'.
//' @return The corresponding results.
//'
//' @examples
//' x <- matrix(rnorm(1e4), 100)
//' y <- matrix(rnorm(1e2), 100)
//' z <- matrix(rnorm(1e4), 100)
//'
//' fMatAdd(x, z) # x + z
//' fMatSubtract(x, z) # x - z
//'
//' fMatElementWiseProduct(x, z) # x * z
//' fMatElementWiseDivide(x, z) # x / z
//'
//' fMatColSum(x) # colSums(x)
//' fMatRowSum(x) # rowSums(x)
//' fMatRowMin(x) # apply(x, 1, min)
//' fMatRowMax(x) # apply(x, 1, max)
//' fMatColMin(x) # apply(x, 2, min)
//' fMatColMax(x) # apply(x, 2, max)
//'
//' fMatSumDiffSquared(x, z) # sum((x-z)^2)
//' @rdname fast_matrix_extra_ops
//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatAdd(SEXP X, SEXP Y) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if ((Rf_nrows(X) != Rf_nrows(Y)) || (Rf_ncols(X) != Rf_ncols(Y))) {
    Rcpp::stop("X and Y must have the same dimensions");
  }

  if ((TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP)) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    Eigen::Map<Eigen::MatrixXi> YMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(Y));
    return Rcpp::wrap(XMti + YMti);
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    return Rcpp::wrap(XMtd + YMtd);
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatSubtract(SEXP X, SEXP Y) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if ((Rf_nrows(X) != Rf_nrows(Y)) || (Rf_ncols(X) != Rf_ncols(Y))) {
    Rcpp::stop("X and Y must have the same dimensions");
  }

  if ((TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP)) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    Eigen::Map<Eigen::MatrixXi> YMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(Y));
    return Rcpp::wrap(XMti - YMti);
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    return Rcpp::wrap(XMtd - YMtd);
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatElementWiseProduct(SEXP X, SEXP Y) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if ((Rf_nrows(X) != Rf_nrows(Y)) || (Rf_ncols(X) != Rf_ncols(Y))) {
    Rcpp::stop("X and Y must have the same dimensions");
  }

  if ((TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP)) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    Eigen::Map<Eigen::MatrixXi> YMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(Y));
    return Rcpp::wrap(XMti.array() * YMti.array());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
    return Rcpp::wrap(XMtd.array() * YMtd.array());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatElementWiseDivide(SEXP X, SEXP Y) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if ((Rf_nrows(X) != Rf_nrows(Y)) || (Rf_ncols(X) != Rf_ncols(Y))) {
    Rcpp::stop("X and Y must have the same dimensions");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
  return Rcpp::wrap(XMtd.array() / YMtd.array());
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatRowSum(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.rowwise().sum());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.rowwise().sum());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatColSum(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.colwise().sum());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.colwise().sum());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatRowMin(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.rowwise().minCoeff());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.rowwise().minCoeff());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatColMin(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.colwise().minCoeff());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.colwise().minCoeff());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatRowMax(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.rowwise().maxCoeff());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.rowwise().maxCoeff());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
SEXP fMatColMax(SEXP X) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP) {
    Eigen::Map<Eigen::MatrixXi> XMti = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(cast_integer(X));
    return Rcpp::wrap(XMti.colwise().maxCoeff());
  } else {
    Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
    return Rcpp::wrap(XMtd.colwise().maxCoeff());
  }
}

//' @name fast_matrix_extra_ops
//' @export
// [[Rcpp::export]]
double fMatSumDiffSquared(SEXP X, SEXP Y) {
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (!(Rf_isMatrix(Y) && (TYPEOF(Y) == REALSXP || TYPEOF(Y) == INTSXP || TYPEOF(Y) == LGLSXP))) {
    Rcpp::stop("'Y' must be a numeric matrix");
  }

  if ((Rf_nrows(X) != Rf_nrows(Y)) || (Rf_ncols(X) != Rf_ncols(Y))) {
    Rcpp::stop("X and Y must have the same dimensions");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  Eigen::Map<Eigen::MatrixXd> YMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(Y));
  return (XMtd-YMtd).array().square().sum();
}
