#include <Rcpp.h>

// [[Rcpp::export]]
SEXP cast_numeric(SEXP input) {
  if (!Rf_isReal(input)) {
    return Rf_coerceVector(input, REALSXP);
  } else {
    return input;
  }
}
