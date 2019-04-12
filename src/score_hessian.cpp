#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
using std::multimap;  
#include <vector>
using std::vector;

#include "bvnorm.h"
#include "utility.h"

// [[Rcpp::plugins(cpp11)]]
// full model 
void score_hessian(arma::vec& score, arma::mat& hessian, // only adds to score and hessian, does not check
                   const vector<double>& rho_v, 
                   const vector<double>& lower_bd, const vector<double>& upper_bd, 
                   const vector<vector<double> >& upper_bdd, 
                   const vector<vector<double> >& lower_bdd, 
                   const multimap<int, vector<int> >& labeled_pairs){
  
  int p_main = upper_bdd[0].size(), p = p_main + rho_v.size(); 
  
  for(int lab = 1; lab != rho_v.size() + 1; ++lab){
    
    double rho = rho_v[lab - 1];
    
    for(auto pos = labeled_pairs.equal_range(lab); pos.first != pos.second; ++pos.first){
      
      vector<int> pair_ind = pos.first->second;
      arma::vec tmp_scr(p); tmp_scr.fill(0); // component score 
      
      for(int pp = 0; pp != p_main; ++pp){
        tmp_scr(pp) = score_main_cpp(pair_ind, pp, rho, lower_bd, upper_bd, upper_bdd, lower_bdd);
      }
      tmp_scr(p_main + lab - 1) = score_rho_cpp(pair_ind, rho, lower_bd, upper_bd, arma::datum::pi);
      
      //update score and hessian
      score += tmp_scr;
      arma::mat tmp_hes(tmp_scr);
      hessian +=  tmp_hes * tmp_hes.t();
    }
  }
 
}


// selected model
void score_hessian(arma::vec& score, arma::mat& hessian, // only adds to score and hessian, does not check
                   const vector<double>& rho_v,  const arma::Col<int>& v_rho, 
                   const int& p_main, const int& p,
                   const vector<double>& lower_bd, const vector<double>& upper_bd, 
                   const vector<vector<double> >& upper_bdd, 
                   const vector<vector<double> >& lower_bdd, 
                   const multimap<int, vector<int> >& labeled_pairs){
  
  for(int lab = 1; lab != rho_v.size() + 1; ++lab){
    
    double rho = rho_v[lab - 1];
    
    for(auto pos = labeled_pairs.equal_range(lab); pos.first != pos.second; ++pos.first){
      
      vector<int> pair_ind = pos.first->second;
      arma::vec tmp_scr(p);  tmp_scr.fill(0);
      
      for(int pp = 0; pp != p_main; ++pp){
        tmp_scr(pp) = score_main_cpp(pair_ind, pp, rho, lower_bd, upper_bd, upper_bdd, lower_bdd);
      }
      
      if(rho){
        arma::Col<int>::const_iterator v_rho_ptr = v_rho.cbegin();
        arma::rowvec::iterator scr_ptr = tmp_scr.begin() + p_main;
        for(int i = 0; i != lab - 1; ++i){
          if(*v_rho_ptr++)
            ++scr_ptr;
        }
        *scr_ptr = score_rho_cpp(pair_ind, rho, lower_bd, upper_bd, arma::datum::pi);
      }
      
      //update score and hessian
      score += tmp_scr;
      arma::mat tmp_hes(tmp_scr);
      hessian += tmp_hes * tmp_hes.t();
    }
  }
}