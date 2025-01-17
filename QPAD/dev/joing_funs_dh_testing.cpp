#include <RcppArmadillo.h>
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppNumerical)]]
using namespace Rcpp;
using namespace RcppArmadillo;
using namespace Numer;


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
class f_d_CPP: public Numer::Func
{
private: 
  double phi_k;
  double tau_k;
  double tmax;
public:
  f_d_CPP( double phi_k_, double tau_k_, double tmax_) : phi_k(phi_k_), tau_k(tau_k_), tmax(tmax_) {}
  
  double operator()(const double& dmax) const
  
  {
    return 2*M_PI*dmax *(1-exp(-phi_k*tmax*exp(-pow(dmax,2)/pow(tau_k,2))));
  }
};

// [[Rcpp::export]]
double nll_fun(NumericVector params, arma::mat X1, arma::mat X2,
             StringVector  tau_params, int nsurvey,
             arma::cube Yarray,
             arma::mat tarray,
             arma::mat rarray,
             NumericVector nrint,
             NumericVector ntint,
             NumericVector max_r,
             NumericVector Ysum,
             NumericVector nlimit
             ){
  NumericVector subset = params[Range(0, tau_params.size()-1)];
  arma::vec sub_v = as<arma::vec>(subset);
  NumericVector subset_phi = params[Range(tau_params.size(), params.size()-1)];
  arma::vec sub_v_phi = as<arma::vec>(subset_phi);


  // arma::vec subset = Rcpp::as<arma::vec>();
  arma::vec tau = exp(X1 * sub_v );
  arma::vec phi = exp(X2 * sub_v_phi );

  NumericVector nll(nsurvey);

  for(int k = 0; k < nsurvey; ++k){
    // std::cout<< "k is " << k  << std::endl;
    const double tau_k = tau[k];
    const double phi_k = phi[k];
    // # Calculate CDF

    arma::mat Y_mat_slice = Yarray.slice(k);//
    // std::cout <<Y_mat.n_rows << std::endl;
    // std::cout << Y_mat.n_cols;
    // // NumericMatrix Y_mat2 = Rcpp::as<Rcpp::NumericMatrix>(Y_mat);
    // arma::mat Ya = Y_mat.submat(0,0, (nrint[k]-1),(ntint[k]-1));
    arma::mat Y_mat = Y_mat_slice(arma::span(0, nrint[k]-1), arma::span(0, ntint[k]-1));
      // Y <- Yarray[k,1:nrint[k],1:ntint[k]]
    NumericMatrix Y = wrap(Y_mat);
    // std::cout<< "Y is " << Y  << std::endl;
    NumericMatrix CDF_binned(nrint[k],ntint[k]);



    for(int j = 0; j < ntint[k]; ++j){
      double tmax = (tarray(k,j)); //max Not sure why R has max of a single value?
      for(int i = 0; i < nrint[k]; ++i){
      // std::cout<< "kji: "<< k <<"\t"<<j<<"\t"<<i << std::endl;
        double upper_r_t = rarray(k,i);
        if(upper_r_t == R_PosInf) upper_r_t = max_r[k];
        const double upper_r = upper_r_t;
        // # Integrate from 1 m from the observer
        f_d_CPP f(phi_k,tau_k,tmax);
        const double lower_limit = 0.01;
        double err_est;
        int err_code;
        
        const double res = integrate(f,lower_limit,
                                    upper_r, err_est, err_code, 500);//,
        // std::cout<< "Result: "<< res << std::endl;
        
         CDF_binned(i,j) = res;
         // CDF_binned(i,j) = 1.2;

                                    // subdiv = 500);

      }
    }
    
    // # Difference across distance bins
    NumericMatrix tmp1(clone(CDF_binned));
        if (tmp1.nrow()>1){
          for (int i=1; i < (nrint[k]);++i){// in 2:nrint[k]){
            for(int j = 0; j < tmp1.ncol(); ++j){
            tmp1(i,j) = CDF_binned(i,j) - CDF_binned(i-1,j);
          }
          }
        }
      // # Difference across time bins
        NumericMatrix p_matrix(clone(tmp1));
        if (p_matrix.ncol()>1){
          for (int j=1; j < (ntint[k]); ++j){// in 2:ntint[k]){
            for(int u=0; u< p_matrix.nrow(); ++u){
              p_matrix(u,j) = tmp1(u,j) - tmp1(u,(j-1));
             
            }
          }
          
        }
        
        double sum_p_matrix =0;
        for (int j=0; j < p_matrix.ncol(); ++j){// in 2:ntint[k]){
          for(int u=0; u< p_matrix.nrow(); ++u){
            sum_p_matrix += p_matrix(u,j);
          }
        }
        
        
        NumericMatrix p_matrix_norm(clone(p_matrix));
        // # Normalize
        for(int i = 0; i < p_matrix.nrow(); ++i ){
          for(int j = 0; j < p_matrix.ncol(); ++j ){
            p_matrix_norm(i,j) = p_matrix(i,j)/sum_p_matrix;
          }
        }
        
        // return(Rcpp::List::create(
        //     Named("CDF_binned")=CDF_binned,
        //     Named("tmp1") = tmp1,
        //     Rcpp::Named("p_matrix") = p_matrix,
        //     Named("normed") = p_matrix_norm));
      // std::cout<<"Y: "<< Y << std::endl;
        // std::cout<< "Ysum: "<< Ysum[k] << std::endl;
      // std::cout<< "Pmatrix: "<<p_matrix_norm << std::endl;
        
        nll[k] = logdmultinomCPP(Y, Ysum[k], p_matrix_norm);
    }// # close loop on k
  // std::cout << "nll"<< nll << std::endl;
      double nll_sum =0;// <- -sum(nll)
  for(int i = 0; i< nsurvey; ++i){
    nll_sum -= nll[i];
  }
  // NumericVector nll_vec = nll_sum;
  
  if ( R_IsNA(nll_sum) ) return(nlimit[2]); 
  // else if (all(is_nan(nll_vec)) ) return(nlimit[2]);
  else if (  !arma::is_finite(nll_sum) ) return(nlimit[2]);
  else return(nll_sum); 
     
    }





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
x = matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 0, 0), ncol=1)
prob = matrix(c(0.2766035, 0.181491 ,0.128348, 0.09668001, 0.07657428, 0.06303872, 0.05344503,0.04634278, 0.04089231, 0.03658438), nrow=1, ncol = 10, byrow = T)
size = 0

logdmultinomCPP(x, size, prob)

# f_d_CPP(1, 1,runif(10, 0, 10),1)
*/
