#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <map>

#include "grid_ring.h"
#include "near_PD.h"

// [[Rcpp::plugins(cpp11)]]
arma::mat cor_mat(const std::vector<double>& rho_v, const std::multimap<int, std::vector<int> >& labeled_pairs, int d){
  arma::mat ret;
            ret.eye(d, d);
  
  for(std::vector<double>::size_type i = 0; i != rho_v.size(); ++i){
    double rho = rho_v[i];
    if(rho){
      for(auto pos = labeled_pairs.equal_range(i + 1); pos.first != pos.second; ++pos.first){
        std::vector<int> pair_ind = pos.first->second;
        ret(pair_ind[0], pair_ind[1]) = rho;
        ret(pair_ind[1], pair_ind[0]) = rho;
      }
    }
  }
  
  nearPD_cpp(ret);
  
  return ret; 
}


