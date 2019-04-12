#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <map>

#include "boot_data.h"
#include "labeled_pairs.h"
#include "make_cor.h"

// [[Rcpp::export]]
Rcpp::List sim_data_cpp(const double y_ini, 
                        const int n, const int K, const int t_size, 
                        const arma::vec& beta, std::vector<double> rho_v){
  
  int p_rho, d = n*n*K; 
  arma::vec y_0(d); y_0.fill(y_ini);
  const std::multimap<int, std::vector<int> > labeled_pairs0 = get_pairs(K, n, p_rho); 
  
  if(p_rho != rho_v.size())
    Rcpp::Rcout << "Incorrect number of correlation parameters!" << std::endl; 
  
  arma::mat corr = cor_mat(rho_v, labeled_pairs0, d); 
  
  return boot_data(y_0, n, K, t_size, beta, corr);
} 
