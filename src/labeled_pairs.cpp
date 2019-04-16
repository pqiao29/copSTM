#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <map>
#include <vector>

#include "grid_ring.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
arma::Mat<int> make_cor_label(int K, int n, int& p_rho, 
                              const int cor_type){
  
  // Block matrix components
  arma::Mat<int> corr_same_loc(K, K, arma::fill::zeros), corr_adj_loc(K, K, arma::fill::zeros); 
  int k = 0; 
  if(cor_type == 1 || cor_type == 3){
    while(k != K){
      corr_adj_loc(k, k) = k + 1;
      ++k;
    } 
  }
  if(cor_type == 2 || cor_type == 3){
    for(int i = 0; i != K - 1; ++i){
      for(int j = i + 1; j != K; ++j){
      corr_same_loc(i, j) = k + 1;
      corr_same_loc(j, i) = k + 1;
      ++k;
     }
   } 
    }
  if(cor_type == 3){ 
    for(int i = 0; i != K - 1; ++i){
      for(int j = i + 1; j != K; ++j){
      corr_adj_loc(i, j) = k + 1;
      corr_adj_loc(j, i) = k + 1;
      ++k;
    }
  } 
 }
  p_rho = k;
  
  std::map<int, std::vector<int> > neighbour = grid_ring(n, 1, false); // neighbour locations
  
  // fill cor_label
  int d = K*n*n;
  arma::Mat<int> ret(d, d, arma::fill::zeros);
  for(int i = 0; i != (n*n); ++i){
    if(cor_type == 2 || cor_type == 3)
       ret.submat(i*K, i*K, (i*K + K - 1), (i*K + K - 1)) = corr_same_loc;
    if(cor_type == 1 || cor_type == 3){
      for(const auto& j : neighbour[i + 1]){
        ret.submat(i*K, (j - 1)*K, (i*K + K - 1), (j*K - 1)) = corr_adj_loc;
        ret.submat((j - 1)*K, i*K, (j*K - 1), (i*K + K - 1)) = corr_adj_loc;
      }
    }
  }
  return ret; 
}

std::multimap<int, std::vector<int> > get_pairs(int K, int n, int& p_rho, 
                                                int rep = 1,  
                                                const int cor_type = 3){
  
  const arma::Mat<int> cor_matrix = make_cor_label(K, n, p_rho, cor_type);
  
  std::multimap<int, std::vector<int> > ret;
    int d = K*n*n;
  for(int i = 0; i != d - 1; ++i){
    for(int j = i + 1; j != d; ++j ){
      int tmp = cor_matrix(i, j);
      if(tmp != 0){
          for(int t = 0; t != rep; ++t){
              ret.insert(make_pair(tmp, std::vector<int> {(t*d + i), (t*d + j)}));
          }
      }
    }
  }
  return ret;
}

