// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
using namespace Numer;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
double logdmultinomCPP(NumericMatrix x,double size, NumericMatrix prob  ) {
  int nrow_M = x.nrow();
  int ncol_M = x.ncol();
  // arma::mat out_mat = x * log(prob ) - lgamma(x + 1);
  double sum_right = 0;
  double left_side = lgamma(size + 1);
  for(int i = 0; i < nrow_M; ++i){
    for(int j = 0; j < ncol_M; ++j){
      double tmp = x(i,j) * log(prob(i,j) ) - lgamma(x(i,j) + 1);
      sum_right += tmp;
    }
  }
    return  left_side + sum_right;
  // lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))

}

// # Calculate CDF and p
// [[Rcpp::export]]


// [[Rcpp::export]]
// List nll_fun(NumericVector params, arma::mat X1, arma::mat X2,
//              StringVector  tau_params, int nsurvey, 
//              arma::cube Yarray ,
//              arma::mat tarray,
//              NumericVector nrint,
//              NumericVector ntint
//              ){
//   NumericVector subset = params[Range(0, tau_params.size()-1)];
//   arma::vec sub_v = as<arma::vec>(subset);
//   NumericVector subset_phi = params[Range(tau_params.size(), params.size()-1)];
//   arma::vec sub_v_phi = as<arma::vec>(subset_phi);
//   
//   
//   // arma::vec subset = Rcpp::as<arma::vec>();
//   arma::vec tau = exp(X1 * sub_v );
//   arma::vec phi = exp(X2 * sub_v_phi );
//   
//   NumericVector nll(nsurvey);
//   
//   for(int k = 0; k < nsurvey; ++k){
//     
//     double tau_k = tau[k];
//     double phi_k = phi[k];
//     // # Calculate CDF
//     
//     arm::Mat Y_mat = Yarray.slice(k);//
//     NumericMatrix Y = Y_mat[Range(0, nrint[k]),Range(0,ntint[k])];
//     NumericMatrix CDF_binned(nrint[k],ntint[k]);
//     
//     NumericVector f_d_CPP( NumericVector dmax){
//       return 2*M_PI*dmax *(1-exp(-phi_k*tmax*exp(-pow(dmax,2)/pow(tau_k,2))));
//     }
//     
//     for(int j = 0; j < ntint[j]; ++j){
//       double tmax = max(tarray[k,j]);
//       for(int i = 0; i < ntint[i]; ++i){
//         double upper_r rarray[k,i];
//         if(upper_r == R_PosInf) upper_r = max_r[k];
//         // # Integrate from 1 m from the observer
//         CDF_binned[i,j] = integrate(f_d_CPP,0.01,
//                                     upper_r,
//                                     subdiv = 500);
//           
//       }
//     }
//     
//     // # Difference across distance bins
//     NumericMatrix tmp1 = CDF_binned;
//         if (nrow(tmp1)>1){
//           for (int i=1; i < nrint[i];++i){// in 2:nrint[k]){
//             tmp1[i,] = CDF_binned[i,] - CDF_binned[i-1,]
//           }
//         }
//         
//     // # Difference across time bins
//         NumericMatrix p_matrix = tmp1;
//         double sum_p_matrix =0;
//         if (ncol(p_matrix)>1){
//           for (int j; j < ntint[k]; ++k){// in 2:ntint[k]){
//             p_matrix[,j] = tmp1[,j] - tmp1[,j-1];
//             for( int r =0; r<nrow(p_matrix); ++r){
//               sum_p_matrix += p_matrix[r,j];
//             }
//           }
//         }
//         
//     // # Normalize
//         p_matrix = p_matrix/sum_p_matrix;
//           
//         nll[k] <- logdmultinomCPP(Y, Ysum[k], p_matrix);
//       }// # close loop on k
//       
//       double nll_sum =0;// <- -sum(nll)
//   for(int i = 0; i< nsurvey; ++i){
//     nll_sum -= nll[i];
//   }
//       // if (nll %in% c(NA, NaN, Inf, -Inf)) return(nlimit[2]) else return(nll);
//       return(nll_sum);
//     }
//       
//       
      
      

  
//   return(
//     Rcpp::List::create(
//       Rcpp::Named("tau") = Rcpp::as<Rcpp::NumericVector>(wrap(tau)),
//       Rcpp::Named("phi") = Rcpp::as<Rcpp::NumericVector>(wrap(phi))
//     ) );
// }



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# x = matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 0, 0), ncol=1)
# prob = matrix(c(0.2766035, 0.181491 ,0.128348, 0.09668001, 0.07657428, 0.06303872, 0.05344503,0.04634278, 0.04089231, 0.03658438), nrow=1, ncol = 10, byrow = T)
# size = 0
# 
# logdmultinomCPP(x, size, prob)

# f_d_CPP(1, 1,runif(10, 0, 10),1)
*/
