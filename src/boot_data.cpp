#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


arma::mat rmvnorm_arma( const arma::vec& mu,
                        const arma::mat& Sigma) {
  unsigned int p = Sigma.n_cols;
  Rcpp::NumericVector draw = Rcpp::rnorm(p);
  arma::mat Z = arma::mat(draw.begin(), 1, p, false, true); 
  arma::mat Y = arma::repmat(mu, 1, 1).t() +  Z * arma::chol(Sigma);
  return Y; 
}

// [[Rcpp::plugins(cpp11)]]
Rcpp::NumericVector rmvpois( const arma::vec& lamvec, const arma::mat& Sigma){
  
  const size_t d = lamvec.size();
  arma::vec tmp_mean(d, arma::fill::zeros); 
  Rcpp::NumericVector XX = Rcpp::wrap(rmvnorm_arma( tmp_mean, Sigma)), ret;
  
  for(size_t i = 0; i != d; ++i){
    double uu = R::pnorm(XX(i), 0, 1, true, false);
    ret.push_back(R::qpois(uu, lamvec(i), true, false));
  }
  return ret;
}


#include<vector>
#include<map>
#include "grid_ring.h"


void get_covariate(arma::mat& x, int& ind_cov_row,
                   const arma::vec& y_prev,   
                   const int& n_lattice, const int& K,
                   std::map<int, std::vector<int> >& nbr){
  // not changing nbr, omit const only because [] is not a const operator
  for(int i = 1; i != n_lattice + 1; ++i){
    int ind_cov_col = 0;
    arma::rowvec tmp_cov(K + 1, arma::fill::zeros); 
    std::vector<int> tmp_nb = nbr[i];
    for(int k = 1; k != K + 1; ++k){
      std::vector<int>::const_iterator iter = tmp_nb.cbegin();
      while(iter != tmp_nb.cend()){
        tmp_cov(k) += log(y_prev((*iter++ - 1)*K + k - 1) + 1);
      }
    } 
    tmp_cov /= tmp_nb.size();
    tmp_cov(0) = 1;
    // fill in covariates
    for(int k = 0; k != K; ++k){ 
      x.submat(ind_cov_row++, ind_cov_col, arma::size(1, K + 1)) = tmp_cov;
      ind_cov_col += K + 1;
    } 
  }
  
} 



Rcpp::List boot_data(const arma::vec& y_0, 
                     const int n, const int K, const int t_size, 
                     const arma::vec& beta, const arma::mat& cor){
  
  /*
   * require: beta.size == K*(K+1), y_0.size() == n*n*K
   */
  
  int n_lattice = n * n, d = n_lattice * K;
  std::map<int, std::vector<int> > nbr = grid_ring(n, 1);
  
  arma::vec y(d * t_size); 
  arma::mat x((d * t_size), (K*(K + 1)), arma::fill::zeros);
  
  // get covariate from y_0 (t = 0)
  int ind_cov_row = 0; 
  get_covariate(x, ind_cov_row, y_0, n_lattice, K, nbr);
  // get t = 1, ..., t_size - 2
  int ind_res = 0; 
  for(int t = 0; t != t_size - 1; ++t){
    arma::vec lmd_t = x.rows(ind_res, (ind_res + d - 1)) * beta;
    arma::vec tmp_y = Rcpp::as<arma::vec>(rmvpois(exp(lmd_t), cor));
    y.subvec(ind_res, (ind_res + d - 1)) = tmp_y;
    ind_res += d;
    get_covariate(x, ind_cov_row, tmp_y, n_lattice, K, nbr);
  }
  arma::vec lmd_t = x.rows(ind_res, (ind_res + d - 1)) * beta;
  y.subvec(ind_res, (ind_res + d - 1)) = Rcpp::as<arma::vec>(rmvpois(exp(lmd_t), cor));

  return Rcpp::List::create(Rcpp::Named("response") = y,
                            Rcpp::Named("covariates") = x);
}


