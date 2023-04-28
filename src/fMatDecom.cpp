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

// [[Rcpp::depends(oneMKL)]]

//' Functions to do the decomposition by leveraging Intel MKL through RcppEigen
//'
//' \describe{
//' \item{\strong{fMatChol}}{This function performs the Cholesky decomposition of the matrix `X`,
//'  i.e., `X = L L^T`, where `L` is a lower triangular matrix with real
//'  and positive diagonal entries, and `L^T` is the transpose of `L`.}
//' \item{\strong{fMatLU}}{This function performs the LU decomposition of the matrix `X`, namely, `X = PLU`,
//'  where `L` is a lower triangular matrix with unit diagonal entries,
//'  `U` is an upper triangular matrix and `P` is a permutation matrix.}
//' \item{\strong{fMatQR}}{This function performs the QR decomposition of the matrix `X`, i.e., `X = QR`,
//'  where `Q` is an orthogonal matrix and `R` is an upper triangular matrix.}
//' \item{\strong{fMatSVD}}{This function performs the singular value decomposition (SVD) of the matrix `X`, namely, `X = U D V^T`,
//'  where `U` and `V` are orthogonal matrices and `D` is a diagonal matrix.}
//' \item{\strong{fMatEigen}}{This function performs the eigenvalue decomposition of the matrix `X`,
//'  i.e., `X = V D V^(-1)`, where 'V' is a matrix whose columns are the eigenvectors of 'X',
//'  and 'D' is a diagonal matrix whose entries are the corresponding eigenvalues of 'X'.}
//' \item{\strong{fMatEigen}}{This function returns a list of objects with two elements:  'values' and 'vectors',
//'  which are respectively the eigenvalues and eigenvectors of 'X'.}
//' }
//'
//' @param X The input matrix.
//' @rdname fast_matrix_decomposition
//' @name fast_matrix_decomposition
//' @examples
//' m <- matrix(c(5,1,1,3),2,2)
//' fMatChol(m)
//' all.equal(fMatChol(m), chol(m)) # It's the same to R
//'
//' X <- matrix(rnorm(9), 3, 3)
//' luRes <- fMatLU(X)
//' solve(luRes$P) %*% luRes$L %*% luRes$U # X = P^(-1) L U
//'
//' qrRes <- fMatQR(X)
//' qrRes$Q %*% qrRes$R # X = Q R
//'
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
//' X <- hilbert(9)[, 1:6]
//' (svdRes <- fMatSVD(X))
//' svdRes$U[ , 1:6] %*% diag(svdRes$d) %*% t(svdRes$V) #  X = U D V'
//' t(svdRes$U[ , 1:6]) %*% X %*% svdRes$V #  D = U' X V
//'
//' X <- hilbert(9)
//' eigenRes <- fMatEigen(X)
//' Re(eigenRes$vectors %*% diag(eigenRes$values) %*% solve(eigenRes$vectors)) # X = V D V^(-1)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatChol(const Eigen::Map<Eigen::MatrixXd> X){
  return X.llt().matrixU();
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatLU(const Eigen::Map<Eigen::MatrixXd> X){
  auto lu = X.partialPivLu();
  auto luMatrix = lu.matrixLU();
  Eigen::MatrixXd L = luMatrix.triangularView<Eigen::StrictlyLower>();
  L.diagonal().setOnes();
  Eigen::MatrixXd U = luMatrix.triangularView<Eigen::Upper>();
  Eigen::MatrixXd P = lu.permutationP();
  return Rcpp::List::create(
    Rcpp::Named("L") = L,
    Rcpp::Named("P") = P,
    Rcpp::Named("U") = U
  );
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatQR(const Eigen::Map<Eigen::MatrixXd> X){
  auto qr = X.householderQr();
  Eigen::MatrixXd Q = qr.householderQ();
  Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q, Rcpp::Named("R") = R
  );
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatSVD(const Eigen::Map<Eigen::MatrixXd> X){
  auto svd = X.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU);
  return Rcpp::List::create(
    Rcpp::Named("d") = svd.singularValues(),
    Rcpp::Named("U") = svd.matrixU(),
    Rcpp::Named("V") = svd.matrixV()
  );
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatEigen(const Eigen::Map<Eigen::MatrixXd> X){
  Eigen::EigenSolver<Eigen::MatrixXd> es(X);
  return Rcpp::List::create(
    Rcpp::Named("values") = es.eigenvalues(), Rcpp::Named("vectors") = es.eigenvectors()
  );
}
