#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
using std::vector;
#include <map>
using std::multimap; using std::map;

#include "grid_ring.h"
#include "near_PD.h"

// [[Rcpp::plugins(cpp11)]]
arma::mat cor_mat(const vector<double>& rho_v, const multimap<int, vector<int> >& labeled_pairs, int d){
  arma::mat ret;
            ret.eye(d, d);
  
  for(vector<double>::size_type i = 0; i != rho_v.size(); ++i){
    double rho = rho_v[i];
    if(rho){
      for(auto pos = labeled_pairs.equal_range(i + 1); pos.first != pos.second; ++pos.first){
        vector<int> pair_ind = pos.first->second;
        ret(pair_ind[0], pair_ind[1]) = rho;
        ret(pair_ind[1], pair_ind[0]) = rho;
      }
    }
  }
  
  nearPD_cpp(ret);
  
  return ret; 
}


/*
#include "labeled_pairs.h"
// [[Rcpp::export]]
arma::mat test(const vector<double>& rho_v, int K, int n){
  //labeled_pairs
  
  int d = K*n*n;
  const multimap<int, vector<int> > labeled_pairs = get_pairs(K, n);
  
  return cor_mat(rho_v, labeled_pairs, d);
}*/