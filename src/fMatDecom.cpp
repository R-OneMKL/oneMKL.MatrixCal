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

//' Functions to do the decomposition by leveraging Intel MKL
//'
//' @param X A matrix to perform decomposition.
//' @rdname fast_matrix_decomposition
//' @name fast_matrix_decomposition
//' @examples
//' m <- matrix(c(5,1,1,3),2,2)
//' fMatChol(m)
//' all.equal(fMatChol(m), chol(m)) # It's the same to R
//' fMatChol(m, FALSE) # lower CHOL matrix
//'
//' X <- matrix(rnorm(9), 3, 3)
//' luRes <- fMatLu(X)
//' solve(luRes$P) %*% luRes$L %*% luRes$U # X = P^(-1) L U
//'
//' qrRes <- fMatQr(X)
//' qrRes$Q %*% qrRes$R # X = Q R
//'
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
//' X <- hilbert(9)[, 1:6]
//' (svdRes <- fMatSvd(X))
//' D <- diag(as.vector(svdRes$d))
//' svdRes$u[ , 1:6] %*% D %*% t(svdRes$v) #  X = U D V'
//' t(svdRes$u[ , 1:6]) %*% X %*% svdRes$v #  D = U' X V
//'
//' X <- hilbert(9)
//' eigenRes <- fMatEigen(X)
//' Re(eigenRes$vectors %*% diag(eigenRes$values) %*% solve(eigenRes$vectors)) # X = V D V^(-1)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd fMatChol(const Eigen::Map<Eigen::MatrixXd> X){
  return X.llt().matrixU();
}

//' @param permutation_matrix Whether the permutation matrix is outputted.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatLu(const Eigen::Map<Eigen::MatrixXd> X){
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
Rcpp::List fMatQr(const Eigen::Map<Eigen::MatrixXd> X){
  auto qr = X.householderQr();
  Eigen::MatrixXd Q = qr.householderQ();
  Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q, Rcpp::Named("R") = R
  );
}

//' @param economical Whether to use economical SVD.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatSvd(const Eigen::Map<Eigen::MatrixXd> X){
  auto svd = X.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU);
  return Rcpp::List::create(
    Rcpp::Named("d") = svd.singularValues(),
    Rcpp::Named("u") = svd.matrixU(),
    Rcpp::Named("v") = svd.matrixV()
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
