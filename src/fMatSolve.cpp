#include <oneMKL.h>
#include <mkl_lapacke.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(oneMKL)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat mkl_real_solve(const arma::mat & x, const arma::mat & b) {
  arma::mat copy_x(x.memptr(), x.n_rows, x.n_cols);
  arma::mat output(b.memptr(), b.n_rows, b.n_cols);
  arma::ivec ipiv(x.n_rows, arma::fill::zeros);
  LAPACKE_dgesv(
    LAPACK_COL_MAJOR, copy_x.n_rows, b.n_cols, copy_x.memptr(),
    copy_x.n_rows, ipiv.memptr(), output.memptr(), copy_x.n_rows
  );
  return output;
}

// [[Rcpp::export]]
arma::cx_mat mkl_cmpl_solve(const arma::cx_mat & x, const arma::cx_mat & b) {
  arma::cx_mat copy_x(x.memptr(), x.n_rows, x.n_cols);
  arma::cx_mat output(b.memptr(), b.n_rows, b.n_cols);
  arma::ivec ipiv(x.n_rows, arma::fill::zeros);
  LAPACKE_zgesv(
    LAPACK_COL_MAJOR, copy_x.n_rows, b.n_cols, copy_x.memptr(),
    copy_x.n_rows, ipiv.memptr(), output.memptr(), copy_x.n_rows
  );
  return output;
}
