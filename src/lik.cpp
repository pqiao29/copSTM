#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
using std::vector;
#include <map>
using std::map;  using std::multimap;

#include "bvnorm.h"
// [[Rcpp::plugins(cpp11)]]
double lik(const vector<double>& rho_v, 
           const multimap<int, vector<int> >& labeled_pairs, 
           const vector<double>& lower_bd, const vector<double>& upper_bd){
  
  double ret = 0;
  
  for(int lab = 1; lab != rho_v.size() + 1; ++lab ){
    for(auto pos = labeled_pairs.equal_range(lab); pos.first != pos.second; ++pos.first){
      
      vector<int> tmp_pair = pos.first->second;
      int i = tmp_pair[0], j = tmp_pair[1];
      
      vector<double> bd1{lower_bd[i], upper_bd[i]}, bd2{lower_bd[j], upper_bd[j]};
      double rho = rho_v[lab - 1];
      ret += log(dbnorm(bd1, bd2, rho));
    }
  }
  return ret;
}
