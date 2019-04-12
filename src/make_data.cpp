#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <vector>

#include "grid_ring.h"
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
Rcpp::List data_cor(const arma::mat& dat, const int n){
  /*
  * Requirement of data columns: Timepoint (0, 1, 2, ...), group (1, 2, 3, ...), tile (1, 2, ... n*n)
  */
  
  std::map<int, std::vector<int> > nbr = grid_ring(n, 1);
  
  int t_size = dat.col(0).max() + 1, K = dat.col(1).max(), n_lattice = n*n;
  
  // count
  arma::cube obs_cnt(K, n_lattice, t_size, arma::fill::zeros);
  for(int i = 0; i != dat.n_rows; ++i){
    ++obs_cnt((dat(i, 1) - 1), (dat(i, 2) - 1), dat(i, 0));
  }
  
  arma::vec response(n_lattice * (t_size - 1) * K); //uninitialized
  arma::mat covariate((n_lattice * (t_size - 1) * K), (K*(K + 1)), arma::fill::zeros);
  
  int ind_cov_row = 0, ind_cov_col, ind_res = 0; 
  for(int t = 1; t != t_size; ++t){
    for(int i = 0; i != n_lattice; ++i){
      //response
      response.subvec(ind_res, (ind_res + K - 1)) = (obs_cnt.slice(t)).col(i); 
      ind_res += K; 
      
      //covariate
      ind_cov_col = 0;
      std::vector<int> tmp_nb = nbr[i + 1]; //neighbour
      arma::rowvec tmp_cov(K + 1, arma::fill::zeros); 
      for(int k = 1; k != K + 1; ++k){
        std::vector<int>::const_iterator iter = tmp_nb.cbegin();
        while(iter != tmp_nb.cend()){
          tmp_cov(k) += log(obs_cnt(k - 1, (*iter++ - 1), t - 1) + 1);
        }
      } // make covariates
      tmp_cov /= tmp_nb.size();
      tmp_cov(0) = 1;
      for(int k = 0; k != K; ++k){
        covariate.submat(ind_cov_row++, ind_cov_col, arma::size(1, K + 1)) = tmp_cov;
        ind_cov_col += K + 1;
      } // fill in covariates
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("covariates") = covariate, 
                            Rcpp::Named("T") = t_size, 
                            Rcpp::Named("K") = K);
  
} 
