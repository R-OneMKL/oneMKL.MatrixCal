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

//' Functions that use oneMKL for fast matrix calculations through RcppEigen
//'
//' \describe{
//' \item{\strong{fMatProd}}{This function returns the multiplication of matrices `X` and `Y`, i.e., `XY`.}
//' \item{\strong{fMatTransProd}}{This function returns the product of the transpose of the matrix `X`
//'  and the matrix `Y`, i.e., `X^T Y`.}
//' \item{\strong{fMatSolve}}{This function returns the solution of a linear system `AX=b`.
//'  If the matrix `X` is symmetric positive definite, Cholesky decomposition
//'  will be used for better computational performance.
//'  If the matrix `X` is invertible, the LU decomposition
//'  will be used for better computational performance.}
//' \item{\strong{fMatInv}}{This function returns the inverse of the matrix `X`, i.e., `X^(-1)`.
//'  If the matrix `X` is symmetric positive definite, Cholesky decomposition
//'  will be used for better computational performance.}
//' \item{\strong{fMatPseudoInv}}{This function returns the pseudo-inverse (also called generalized inverse or g-inverse) of the matrix 'X'.}
//' \item{\strong{fMatLeastSquare}}{This function returns the solution of least square through QR decomposition.
//' If `stable` = TRUE, a more stable but slower algorithm will be applied.
//' If `is_X_full_rank` = FALSE, the pseudo-inverse will be applied.}
//' \item{\strong{fMatAdd}}{This function returns the sum of matrices `X` and `Y`, i.e., `X + Y`}
//' \item{\strong{fMatSubtract}}{This function returns the result of the matrix `X` minus the matrix `Y`, namely, `X - Y`.}
//' \item{\strong{fMatDet}}{This function returns the determinant of the matrix `X`.}
//' \item{\strong{fMatRank}}{This function returns the rank of the matrix `X`.}
//' \item{\strong{fMatRowSum}}{This function returns the sum of each row.}
//' \item{\strong{fMatColSum}}{This function returns the sum of each column.}
//' \item{\strong{fMatRowMin}}{This function returns the minimum of each row.}
//' \item{\strong{fMatRowMax}}{This function returns the maximum of each row.}
//' \item{\strong{fMatColMin}}{This function returns the minimum of each column.}
//' \item{\strong{fMatColMax}}{This function returns the maximum of each column.}
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
//' fMatAdd(x, z) # x + z
//' fMatSubtract(x, z) # x - z
//'
//' fMatDet(x)
//' fMatRank(x)
//'
//' fMatColSum(x) # colSums(x)
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
    Rcpp::NumericMatrix Y,
    bool is_X_symmetric = false
) {
  if (Y.ncol() == 1) {
    // to have the best performance
    Eigen::Map<Eigen::VectorXd> b(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(Y));
    if (is_X_symmetric) {
      return X.selfadjointView<Eigen::Lower>() * b;
    } else {
      return X * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Y));
    if (is_X_symmetric) {
      return X.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return X * Z;
    }
  }
}

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatTransProd(
    const Eigen::Map<Eigen::MatrixXd> X,
    Rcpp::NumericMatrix Y,
    bool is_X_symmetric = false
) {
  if (Y.cols() == 1) {
    // to have best performance
    Eigen::Map<Eigen::VectorXd> b(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(Y));
    if (is_X_symmetric) {
      return X.selfadjointView<Eigen::Lower>() * b;
    } else {
      return X.transpose() * b;
    }
  } else {
    Eigen::Map<Eigen::MatrixXd> Z(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Y));
    if (is_X_symmetric) {
      return X.selfadjointView<Eigen::Lower>() * Z;
    } else {
      return X.transpose() * Z;
    }
  }
}

//' @param is_invertible A logical variable indicating whether the input matrix `X` is invertible.
//'  Better computational performance is expected if the matrix is invertible.
//' @param is_sym_pd A logical variable indicating whether the input matrix `X` is symmetric positive definitive.
//'  Better computational performance is expected if the matrix is symmetric positive definitive.
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
    return X.householderQr().solve(Y);
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
Eigen::MatrixXd fMatPseudoInv(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.completeOrthogonalDecomposition().pseudoInverse();
}

//' @param stable A logical variable indicating whether to use a more stable but slower algorithm.
//' @param is_X_full_rank A logical variable indicating whether the input matrix 'X' is full-rank.
//'  If false, the pseudo inverse will be applied.
//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatLeastSquare(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::VectorXd> Y,
    bool stable = true,
    bool is_X_full_rank = true
) {
  if (is_X_full_rank) {
    if (stable) {
      return X.householderQr().solve(Y);
    } else {
      return X.colPivHouseholderQr().solve(Y);
    }
  } else {
    return X.completeOrthogonalDecomposition().solve(Y);
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

//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatRank(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.colPivHouseholderQr().rank();
}


//' @name fast_matrix_ops
//' @export
// [[Rcpp::export]]
double fMatRCond(const Eigen::Map<Eigen::MatrixXd> X) {
  return X.partialPivLu().rcond();
}
