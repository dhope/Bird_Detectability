#include <Rcpp.h>

// [[Rcpp::export]]
bool isReallyNA(Rcpp::NumericVector val) {
  return Rcpp::all(Rcpp::is_na(val));
}

/*** R
isReallyNA(1)
isReallyNA(1L)
isReallyNA(NA)
isReallyNA(NaN)
isReallyNA(Inf)
*/