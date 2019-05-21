#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "boot_data.h"
#include "labeled_pairs.h"
#include "make_cor.h"

// [[Rcpp::export]]
Rcpp::List sim_data_cpp(const int n, const int K, 
                        bool temporal, const int t_size, 
                        const arma::vec& beta, std::vector<double> rho_v,
                        const int cor_type, const double y_ini, 
                        int marginal, double dispersion){
  
  int p_rho, d = n*n*K; 
  arma::vec y_0(d); y_0.fill(y_ini);
  const std::multimap<int, std::vector<int> > labeled_pairs0 = get_pairs(K, n, p_rho, 1, cor_type); 
  
  if(p_rho != rho_v.size())
    throw Rcpp::exception("Incorrect number of correlation parameters!", false);
  
  arma::mat corr = cor_mat(rho_v, labeled_pairs0, d); 
  
  if(temporal){
    return boot_data(y_0, n, K, t_size, beta, corr);
  }else{
    return gen_data(t_size, d, beta, marginal, dispersion, corr);
  }
} 


