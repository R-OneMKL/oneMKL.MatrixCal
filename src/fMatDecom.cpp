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
#include <string>

// [[Rcpp::depends(oneMKL)]]

//' Functions to do the decomposition by leveraging Intel MKL
//'
//' @param x A matrix to perform decomposition.
//' @param upper A Boolean value to indicate the output matrix is a upper matrix. False will return a lower matrix.
//' @rdname fast_matrix_decomposition
//' @name fast_matrix_decomposition
//' @examples
//' m <- matrix(c(5,1,1,3),2,2)
//' fMatChol(m)
//' all.equal(fMatChol(m), chol(m)) # It's the same to R
//' fMatChol(m, FALSE) # lower CHOL matrix
//'
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
//' X <- hilbert(9)[, 1:6]
//' (s <- fMatSvd(X))
//' D <- diag(as.vector(s$d))
//' s$u[ , 1:6] %*% D %*% t(s$v) #  X = U D V'
//' t(s$u[ , 1:6]) %*% X %*% s$v #  D = U' X V
//'
//' fMatEigen(cbind(c(1,-1), c(-1,1)), TRUE)
//' fMatEigen(cbind(c(1,-1), c(-1,1)), FALSE) # Same, but different datatype
//'
//' X <- matrix(rnorm(9), 3, 3)
//' res <- fMatLu(X, TRUE)
//' # Note that L is generally not lower-triangular when permutation_matrix = FALSE
//' res$P %*% res$L %*% res$U # X = P' L U
//'
//' schurRes <- fMatSchur(X)
//' # Note that Schur decomposition is not unique in general.
//' schurRes$U %*% schurRes$S %*% t(schurRes$U) # X = U S U'
//'
//' qrRes <- fMatQr(X)
//' qrRes$Q %*% qrRes$R # X = Q R
//' @export
// [[Rcpp::export]]
arma::mat fMatChol(const arma::mat & x, bool upper = true){
  std::string layout = upper?"upper":"lower";
  return arma::chol(x, layout.c_str());
}

//' @param economical Whether to use economical SVD.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatSvd(const arma::mat & x, bool economical = false){
  arma::vec d;
  arma::mat u, v;
  if (economical) {
    arma::svd_econ(u, d, v, x);
  } else {
    arma::svd(u, d, v, x);
  }
  return Rcpp::List::create(
    Rcpp::Named("d") = d, Rcpp::Named("u") = u, Rcpp::Named("v") = v
  );
}

//' @param is_symmetric Whether the matrix is symmetric.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatEigen(const arma::mat & x, bool is_symmetric = false){
  if (is_symmetric) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, x);
    return Rcpp::List::create(
      Rcpp::Named("values") = eigval, Rcpp::Named("vectors") = eigvec
    );
  } else {
    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, x);
    return Rcpp::List::create(
      Rcpp::Named("values") = eigval, Rcpp::Named("vectors") = eigvec
    );
  }
}

//' @param permutation_matrix Whether the permutation matrix is outputted.
//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatLu(const arma::mat & x, bool permutation_matrix = false){
  if (permutation_matrix) {
    arma::mat L, U, P;
    arma::lu(L, U, P, x);
    return Rcpp::List::create(
      Rcpp::Named("L") = L, Rcpp::Named("P") = P, Rcpp::Named("U") = U
    );
  } else {
    arma::mat L, U;
    arma::lu(L, U, x);
    return Rcpp::List::create(
      Rcpp::Named("L") = L, Rcpp::Named("U") = U
    );
  }
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatSchur(const arma::mat & x){
  arma::mat U, S;
  arma::schur(U, S, x);
  return Rcpp::List::create(
    Rcpp::Named("U") = U, Rcpp::Named("S") = S
  );
}

//' @name fast_matrix_decomposition
//' @export
// [[Rcpp::export]]
Rcpp::List fMatQr(const arma::mat & x, bool permutation_matrix = false, bool economical = false){
  if (permutation_matrix) {
    arma::mat Q, R;
    arma::umat P;
    arma::qr(Q, R, P, x);
    return Rcpp::List::create(
      Rcpp::Named("Q") = Q, Rcpp::Named("P") = P, Rcpp::Named("R") = R
    );
  } else {
    arma::mat Q, R;
    if (economical) {
      arma::qr_econ(Q, R, x);
    } else {
      arma::qr(Q, R, x);
    }
    return Rcpp::List::create(
      Rcpp::Named("Q") = Q, Rcpp::Named("R") = R
    );
  }
}
