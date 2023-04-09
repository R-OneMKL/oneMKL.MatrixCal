#include <oneMKL.h>
#include <mkl_lapacke.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(oneMKL, RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat mkl_real_solve(const arma::mat & a, const arma::mat & b, double tol) {
  if (!a.is_square()) {
    Rcpp::stop("a (%d x %d) must be square.", a.n_rows, a.n_cols);
  }

  if ((a.n_rows == 0) || (a.n_cols == 0)) {
    Rcpp::stop("Any dimention of a must be greater than 0.");
  }

  if (b.n_rows != a.n_cols) {
    Rcpp::stop("b (%d x %d) must be compatible with a (%d x %d)", b.n_rows, b.n_cols, a.n_rows, a.n_cols);
  }

  if ((b.n_rows == 0) || (b.n_cols == 0)) {
    Rcpp::stop("Any dimention of b must be greater than 0.");
  }

  // initialize ipiv
  arma::ivec ipiv(a.n_rows, arma::fill::zeros);

  // make a copy of a
  arma::mat a2(a.memptr(), a.n_rows, a.n_cols);

  // make a copy of b as output
  arma::mat output(b.memptr(), b.n_rows, b.n_cols);

  // employ LAPACKE_dgesv to solve
  int info = LAPACKE_dgesv(
    LAPACK_COL_MAJOR, a2.n_rows, b.n_cols, a2.memptr(),
    a2.n_rows, ipiv.memptr(), output.memptr(), a2.n_rows
  );

  if (info < 0) {
    Rcpp::stop("argument %d of Lapacke function %s had invalid value", -info, "LAPACKE_dgesv");
  }

  if (info > 0) {
    Rcpp::stop("Lapacke function %s: system is exactly singular: U[%d, %d] = 0", "LAPACKE_dgesv", info, info);
  }

  if (tol > 0.0) {
    char one = '1';
    double anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, one, a.n_rows, a.n_cols, a.memptr(), a.n_rows);
    double rcond = 0.0;
    int info = LAPACKE_dgecon(LAPACK_COL_MAJOR, one, a.n_rows, a.memptr(), a.n_rows, anorm, &rcond);
    if (rcond < tol) {
      Rcpp::stop("system is computationally singular: reciprocal condition number = %g", rcond);
    }
  }

  return output;
}

// [[Rcpp::export]]
arma::cx_mat mkl_cmpl_solve(const arma::cx_mat & a, const arma::cx_mat & b) {
  if (!a.is_square()) {
    Rcpp::stop("a (%d x %d) must be square.", a.n_rows, a.n_cols);
  }

  if ((a.n_rows == 0) || (a.n_cols == 0)) {
    Rcpp::stop("Any dimention of a must be greater than 0.");
  }

  if (b.n_rows != a.n_cols) {
    Rcpp::stop("b (%d x %d) must be compatible with a (%d x %d)", b.n_rows, b.n_cols, a.n_rows, a.n_cols);
  }

  if ((b.n_rows == 0) || (b.n_cols == 0)) {
    Rcpp::stop("Any dimention of b must be greater than 0.");
  }


  // initialize ipiv
  arma::ivec ipiv(a.n_rows, arma::fill::zeros);

  // make a copy of a
  arma::cx_mat a2(a.memptr(), a.n_rows, a.n_cols);

  // make a copy of b as output
  arma::cx_mat output(b.memptr(), b.n_rows, b.n_cols);

  // employ LAPACKE_zgesv to solve
  int info = LAPACKE_zgesv(
    LAPACK_COL_MAJOR, a2.n_rows, b.n_cols, a2.memptr(),
    a2.n_rows, ipiv.memptr(), output.memptr(), a2.n_rows
  );

  if (info < 0) {
    Rcpp::stop("argument %d of Lapacke function %s had invalid value", -info, "LAPACKE_zgesv");
  }

  if (info > 0) {
    Rcpp::stop("Lapacke function %s: system is exactly singular: U[%d, %d] = 0", "LAPACKE_zgesv", info, info);
  }

  return output;
}
