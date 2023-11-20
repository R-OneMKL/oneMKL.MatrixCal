// Copyright (C) 2022             Ching-Chuan Chen
//
// This file is part of oneMKL.
//
// oneMKL.MatrixCal is free software: you can redistribute it and/or
// modify itunder the terms of the GNU General Public License as
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

//' Functions to do the decomposition by leveraging Intel MKL through RcppEigen
//'
//' \describe{
//' \item{\strong{fMatChol}}{This function performs the Cholesky decomposition of the matrix `X`,
//'  i.e., `X = L L^T`, where `L` is a lower triangular matrix with real
//'  and positive diagonal entries, and `L^T` is the transpose of `L`.}
//' \item{\strong{fMatLU}}{This function performs the LU decomposition of the matrix `X`, namely, `PX = LU`,
//'  where `L` is a lower triangular matrix with unit diagonal entries,
//'  `U` is an upper triangular matrix and `P` is a permutation matrix.}
//' \item{\strong{fMatQR}}{This function performs the QR decomposition of the matrix `X`, i.e., `X = QR`,
//'  where `Q` is an orthogonal matrix and `R` is an upper triangular matrix.
//'  If `with_permutation_matrix` = TRUE, it returns a permutation matrix `P` as well.
//'  The decomposition will be `XP = QR` which the algorithm is more stable than the one without a permutation matrix.}
//' \item{\strong{fMatEigen}}{This function performs the eigenvalue decomposition of the matrix `X`,
//'  i.e., `X = V D V^(-1)`, where 'V' is a matrix whose columns are the eigenvectors of 'X',
//'  and 'D' is a diagonal matrix whose entries are the corresponding eigenvalues of 'X'.
//'  This function returns a list of objects with two elements: `values` and `vectors`,
//'  which are respectively the eigenvalues and eigenvectors of `X`.
//'  If `is_X_symmetric` = TRUE, only the real part will be returned.}
//' \item{\strong{fMatSVD}}{This function performs the singular value decomposition (SVD) of the matrix `X`, namely, `X = U D V^T`,
//'  where `U` and `V` are orthogonal matrices and `D` is a diagonal matrix.}
//' }
//'
//' @param X The input matrix.
//' @rdname fast_matrix_decomposition
//' @name fast_matrix_decomposition
//' @examples
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
//' # Cholesky decomposition
//' X <- hilbert(9)
//' fMatChol(X)
//' all.equal(fMatChol(X), chol(X)) # It's the same to R
//'
//' # LU Decomposition
//' luRes <- fMatLU(X)
//' fMatInv(luRes$P) %*% luRes$L %*% luRes$U # X = P^(-1) L U
//'
//' # QR Decomposition
//' qrRes <- fMatQR(X)
//' qrRes$Q %*% qrRes$R # X = Q R
//'
//' # QR Decomposition with Permutation matrix
//' qrRes <- fMatQR(X, TRUE)
//' qrRes$Q %*% qrRes$R %*% fMatInv(qrRes$P) # X = Q R P^{-1}
//'
//' # Eigen Decomposition
//' eigenRes <- fMatEigen(X)
//' Re(eigenRes$vectors %*% diag(eigenRes$values) %*% solve(eigenRes$vectors)) # X = V D V^(-1)
//'
//' # Eigen Decomposition for the symmetric matrix
//' eigenRes2 <- fMatEigen(X, TRUE)
//' eigenRes2$vectors %*% diag(eigenRes2$values) %*% fMatInv(eigenRes2$vectors)
//'
//' # SVD Decomposition
//' Z <- hilbert(9)[, 1:6]
//' (svdRes <- fMatSVD(Z))
//' svdRes$U[ , 1:6] %*% diag(svdRes$d) %*% t(svdRes$V) #  Z = U D V'
//' t(svdRes$U[ , 1:6]) %*% Z %*% svdRes$V #  D = U' Z V
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatChol(SEXP X){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }

  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  return XMtd.llt().matrixU();
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatLU(SEXP X){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  auto lu = XMtd.partialPivLu();
  auto luMatrix = lu.matrixLU();
  Eigen::MatrixXd L = luMatrix.triangularView<Eigen::StrictlyLower>();
  L.diagonal().setOnes();
  Eigen::MatrixXd U = luMatrix.triangularView<Eigen::Upper>();
  Eigen::MatrixXd P = lu.permutationP();
  return Rcpp::List::create(
    Rcpp::Named("P") = P,
    Rcpp::Named("L") = L,
    Rcpp::Named("U") = U
  );
}

//' @param with_permutation_matrix A logical variable indicating whether
//' the QR decomposition is performed with a permutation matrix `P` returned.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatQR(SEXP X, bool with_permutation_matrix = false){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (with_permutation_matrix) {
    auto pqr = XMtd.colPivHouseholderQr();
    Eigen::MatrixXd Q = pqr.householderQ();
    Eigen::MatrixXd R = pqr.matrixQR().triangularView<Eigen::Upper>();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType PMat(pqr.colsPermutation());
    Eigen::MatrixXd P = PMat * Eigen::MatrixXd::Identity(PMat.cols(), PMat.cols());
    return Rcpp::List::create(
      Rcpp::Named("P") = P, Rcpp::Named("Q") = Q, Rcpp::Named("R") = R
    );
  } else {
    auto qr = XMtd.householderQr();
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    return Rcpp::List::create(
      Rcpp::Named("Q") = Q, Rcpp::Named("R") = R
    );
  }
}

//' @param is_X_symmetric A logical variable indicating whether the input matrix `X` is symmetric.
//'  If true, `fMatEigen` returns a real matrix. Otherwise, it returns a complex matrix.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatEigen(SEXP X, bool is_X_symmetric = false){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }

  if (Rf_nrows(X) != Rf_ncols(X)) {
    Rcpp::stop("X must be a square matrix");
  }


  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  if (is_X_symmetric) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(XMtd);
    return Rcpp::List::create(
      Rcpp::Named("values") = es.eigenvalues(), Rcpp::Named("vectors") = es.eigenvectors()
    );
  } else {
    Eigen::EigenSolver<Eigen::MatrixXd> es(XMtd);
    return Rcpp::List::create(
      Rcpp::Named("values") = es.eigenvalues(), Rcpp::Named("vectors") = es.eigenvectors()
    );
  }
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatSVD(SEXP X){
  if (!(Rf_isMatrix(X) && (TYPEOF(X) == REALSXP || TYPEOF(X) == INTSXP || TYPEOF(X) == LGLSXP))) {
    Rcpp::stop("'X' must be a numeric matrix");
  }
  Eigen::Map<Eigen::MatrixXd> XMtd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(cast_numeric(X));
  auto svd = XMtd.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU);
  return Rcpp::List::create(
    Rcpp::Named("d") = svd.singularValues(),
    Rcpp::Named("U") = svd.matrixU(),
    Rcpp::Named("V") = svd.matrixV()
  );
}
