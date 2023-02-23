#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int runit(int nsurveys, arma::cube Yarray, Rcpp::NumericVector nrint, Rcpp::NumericVector ntint){
  for(int k = 0; k < nsurveys; ++k){
    // arma::mat Ymat = Yarray.slice(k);
    // arma::mat tmpMat = Ymat.rows(0, nrint[k]-1);
    // arma::mat Y = tmpMat.cols(0, ntint[k]-1 );
    arma::mat Y(nrint[k],ntint[k] );
    for(int i = 0; i < nrint[k]; ++i){
      for(int j = 0; j < ntint[k]; ++j){
        int z = Yarray(i,j,k);
        Y(i,j) = z;
      }
    }
    }
  return(1);
}

// [[Rcpp::export]]
int runit2(int nsurveys, arma::cube Yarray, Rcpp::NumericVector nrint, Rcpp::NumericVector ntint){
  for(int k = 0; k < nsurveys; ++k){
    arma::mat Y(nrint[k],ntint[k] );
    for(int i = 0; i < nrint[k]; ++i){
      for(int j = 0; j < ntint[k]; ++j){
        Y(i,j) =  Yarray(i,j,k);;
      }
    }
  }
  return(1);
}


// [[Rcpp::export]]
int runit_slice(int nsurveys, arma::cube Yarray, Rcpp::NumericVector nrint, Rcpp::NumericVector ntint){
  for(int k = 0; k < nsurveys; ++k){
    arma::mat Ymat = Yarray.slice(k);
    arma::mat tmpMat = Ymat.rows(0, nrint[k]-1);
    arma::mat Y = tmpMat.cols(0, ntint[k]-1 );
    
  }
  return(1);
}


// Rcpp::NumericMatrix Y(nrint[k],ntint[k] );
// for(int i = 0; i < nrint[k]; ++i){
//   for(int j = 0; j < ntint[k]; ++j){
//     Y(i,j) = Yarray(i,j,k);
//   }
// }

// arma::mat Ymatcl = Ymat(arma::span(0, nrint[k]-1),
//                         arma::span(0, ntint[k]-1));
// Rcpp::NumericMatrix Y = Rcpp::wrap(Ymatcl);

/*** R
load("E:/GitHUB_Clones/Bird_Detectability/QPAD/dev/BreakingBad.RData")
Yarray_x = aperm(Yarray, c(2,3,1))




runit(nsurvey, Yarray_x, nrint, ntint)
bench::mark(c())
out <- purrr::map_dbl(1:20000,
                      ~{print(.x);runit(nsurvey, Yarray_x, nrint, ntint)})

library("microbenchmark")
port_timings <- microbenchmark(r_orig    = runit(nsurvey, Yarray_x, nrint, ntint),
                               r_orig2    = runit2(nsurvey, Yarray_x, nrint, ntint),
                               r_cached  = runit_slice(nsurvey, Yarray_x, nrint, ntint),
                               times = 100000)




*/
