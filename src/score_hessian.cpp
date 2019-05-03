#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bvnorm.h"
#include "utility.h"

// [[Rcpp::plugins(cpp11)]]
// full model 
void score_hessian(arma::vec& score, arma::mat& hessian, // only adds to score and hessian, does not check
                   const std::vector<double>& rho_v, 
                   const std::vector<double>& lower_bd, const std::vector<double>& upper_bd, 
                   const std::vector<std::vector<double> >& upper_bdd, 
                   const std::vector<std::vector<double> >& lower_bdd, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs){
  
  int p_main = upper_bdd[0].size(), p = p_main + rho_v.size(); 
  
  for(int lab = 1; lab != rho_v.size() + 1; ++lab){
    
    double rho = rho_v[lab - 1];
    
    for(auto pos = labeled_pairs.equal_range(lab); pos.first != pos.second; ++pos.first){
      
      std::vector<int> pair_ind = pos.first->second;
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
                   const std::vector<double>& rho_v,  const arma::Col<int>& v_rho, 
                   const int& p_main, const int& p,
                   const std::vector<double>& lower_bd, const std::vector<double>& upper_bd, 
                   const std::vector<std::vector<double> >& upper_bdd, 
                   const std::vector<std::vector<double> >& lower_bdd, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs){
  
  for(int lab = 1; lab != rho_v.size() + 1; ++lab){
    
    double rho = rho_v[lab - 1];
    
    for(auto pos = labeled_pairs.equal_range(lab); pos.first != pos.second; ++pos.first){
      
      std::vector<int> pair_ind = pos.first->second;
      arma::vec tmp_scr(p);  tmp_scr.fill(0);
      
      for(int pp = 0; pp != p_main; ++pp){
        tmp_scr(pp) = score_main_cpp(pair_ind, pp, rho, lower_bd, upper_bd, upper_bdd, lower_bdd);
      }
      
      if(v_rho[lab - 1]){
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

#include "bounds.h"
#include "bounds_deriv.h"
// for obtain standard error in non-temporal setting
// full model
double get_std_err(const arma::mat& xx, const arma::vec& y, 
                   const arma::vec& beta, const std::vector<double>& rho_v, 
                   int t_size, int d, arma::vec& se, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs0){
  
  int p_main = beta.size(), p_rho = rho_v.size(), p = p_main + p_rho;
  arma::vec score(p, arma::fill::zeros); 
  arma::mat hessian(p, p, arma::fill::zeros); 
  arma::mat score_record(p, t_size);
  
  // Bounds for integrals in Gaussian cdf
  const arma::vec y_minus1 = y - 1;
  std::vector<double> lower_bd = bound(xx, y_minus1, beta), 
                      upper_bd = bound(xx, y, beta);
  
  std::vector<std::vector<double> > upper_bdd = get_bd_deriv_arma(xx, y, true, beta),
                                    lower_bdd = get_bd_deriv_arma(xx, y, false, beta);
  
  // score, score_record, hessian 
  for(int t = 0; t != t_size; ++t){
    
    arma::vec score_t(p, arma::fill::zeros);
    
    for(int lab = 1; lab != rho_v.size() + 1; ++lab){
      
      double rho = rho_v[lab - 1];
      
      for(auto pos = labeled_pairs0.equal_range(lab); pos.first != pos.second; ++pos.first){
        
        std::vector<int> pair_ind = pos.first->second;
                         pair_ind[0] += t*d; pair_ind[1] += t*d;
        arma::vec tmp_scr(p); tmp_scr.fill(0); // component score 
        
        for(int pp = 0; pp != p_main; ++pp){
          tmp_scr(pp) = score_main_cpp(pair_ind, pp, rho, lower_bd, upper_bd, upper_bdd, lower_bdd);
        }
        tmp_scr(p_main + lab - 1) = score_rho_cpp(pair_ind, rho, lower_bd, upper_bd, arma::datum::pi);
        
        //update score and hessian
        score += tmp_scr;  
        score_t += tmp_scr;
        arma::mat tmp_hes(tmp_scr);
        hessian +=  tmp_hes * tmp_hes.t();
      }
    }
    
    score_record.col(t) = score_t;
  }
  
  // J
  score_record.each_col() -= (score/t_size);
  arma::mat J = (score_record * score_record.t());
  
  // standard error
  arma::mat tmp_mat = arma::solve(hessian, J);
  arma::mat se_mat = tmp_mat * arma::inv(hessian);
  se = arma::diagvec(se_mat) * t_size;

  return arma::trace(tmp_mat);
}

// selected model
double get_std_err(const arma::mat& xx, const arma::vec& y, 
                   const arma::vec& beta, std::vector<double> rho_v, const int t_size, 
                   const arma::Col<int>& v_main, const arma::Col<int>& v_rho, 
                   int d, arma::vec& se, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs0){
  // set up
  int p_main = sum(v_main), p = p_main + sum(v_rho);
  arma::vec beta_sub = beta % v_main;
  arma::vec rho_zeros = arma::vec(rho_v) % v_rho;
  rho_v = arma::conv_to<std::vector<double> >::from(rho_zeros);
  
  // bounds 
  const arma::vec y_minus1 = y - 1;
  std::vector<double> lower_bd = bound(xx, y_minus1, beta_sub), 
                      upper_bd = bound(xx, y, beta_sub);
  
  std::vector<std::vector<double> > upper_bdd = get_bd_deriv_arma_sub(xx, y, true, beta_sub, v_main),
                                    lower_bdd = get_bd_deriv_arma_sub(xx, y, false, beta_sub, v_main);
  
  // score, hessian, score_record 
  arma::vec score(p, arma::fill::zeros); 
  arma::mat hessian(p, p, arma::fill::zeros); 
  arma::mat score_record(p, t_size);
  
  for(int t = 0; t != t_size; ++t){
    
    arma::vec score_t(p, arma::fill::zeros);
    
    for(int lab = 1; lab != rho_v.size() + 1; ++lab){
      
      double rho = rho_v[lab - 1];
      
      for(auto pos = labeled_pairs0.equal_range(lab); pos.first != pos.second; ++pos.first){
        
        std::vector<int> pair_ind = pos.first->second;
                         pair_ind[0] += t*d; pair_ind[1] += t*d;
        arma::vec tmp_scr(p);  tmp_scr.fill(0);
        
        for(int pp = 0; pp != p_main; ++pp){
          tmp_scr(pp) = score_main_cpp(pair_ind, pp, rho, lower_bd, upper_bd, upper_bdd, lower_bdd);
        }
        
        if(v_rho[lab - 1]){
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
        score_t += tmp_scr;
        arma::mat tmp_hes(tmp_scr);
        hessian += tmp_hes * tmp_hes.t();
      }
    }
    
    score_record.col(t) = score_t;
  }
  
  // J
  score_record.each_col() -= (score/t_size);
  arma::mat J = (score_record * score_record.t());
  
  // standard error
  arma::mat tmp_mat = arma::solve(hessian, J);
  arma::mat se_mat = tmp_mat * arma::inv(hessian);
  se = arma::diagvec(se_mat) * t_size;
  
  return arma::trace(tmp_mat);
}

