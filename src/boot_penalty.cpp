#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "boot_data.h"
#include "bounds.h"
#include "bounds_deriv.h"
#include "score_hessian.h"

#include <map>
using std::multimap;
#include <vector>
using std::vector;


//withstandard error
double boot_CLIC_penalty(  const arma::vec& y_0,
                           const int n, const int K, const int t_size, 
                           const arma::vec& beta, arma::vec& se, 
                           const arma::mat& cor,   // cor required to be valid correlation matrix! 
                           const vector<double>& rho_v,
                           const int& p_main, const int& p,
                           int marginal, double dispersion,
                           const multimap<int, vector<int> >& labeled_pairs,
                           int B, bool Message_prog){
  
  if(Message_prog) Rcpp::Rcout << "Bootstrapping for full model: ";
  
  arma::mat H(p, p, arma::fill::zeros), J(p, p, arma::fill::zeros); 
  int prt = B/20 + 1;
  for(int b = 0; b != B; ++b){
    
    if(Message_prog && (b % prt == 0)) Rcpp::Rcout << "==";
    
    Rcpp::List data = boot_data(y_0, n, K, t_size, beta, cor, marginal, dispersion);
    const arma::vec &y = data["response"];
    const arma::mat &xx = data["covariates"];
    
    const arma::vec y_minus1 = y - 1;
    vector<double> lower_bd = bound(xx, y_minus1, beta, marginal, dispersion), 
                   upper_bd = bound(xx, y, beta, marginal, dispersion);
    
    vector<vector<double> > upper_bdd = get_bd_deriv_arma(xx, y, true, beta, marginal, dispersion),
                            lower_bdd = get_bd_deriv_arma(xx, y, false, beta, marginal, dispersion);
    
    arma::vec tmp_score(p, arma::fill::zeros);
    score_hessian(tmp_score, H, rho_v, lower_bd, upper_bd, 
                  upper_bdd, lower_bdd, labeled_pairs, marginal);
    J += tmp_score * tmp_score.t();
  }
  H /= B;  J /= B; 
  
  if(Message_prog) Rcpp::Rcout << std::endl;
  
  arma::mat tmp_mat = arma::solve(H, J);
  
  arma::mat se_mat = tmp_mat * arma::inv(H);
  se = arma::diagvec(se_mat);
  
  return arma::trace(tmp_mat);
}


// with standard error 
double boot_CLIC_penalty_sub(  const arma::vec& y_0,
                           const int n, const int K, const int t_size, 
                           const arma::vec& beta, arma::vec& se, 
                           const arma::mat& cor, 
                           const vector<double>& rho_v,
                           const arma::Col<int>& v_main, const arma::Col<int>& v_rho,
                           const int& p_main_sub, const int& p_sub,
                           int marginal, double dispersion,
                           const multimap<int, vector<int> >& labeled_pairs, 
                           int B, bool Message_prog){
  
  arma::mat H(p_sub, p_sub, arma::fill::zeros), J(p_sub, p_sub, arma::fill::zeros);
  arma::vec beta_sub = beta % v_main;
  
  int prt = B/20 + 1;
  for(int b = 0; b != B; ++b){
    
    if(Message_prog && (b % prt == 0)) Rcpp::Rcout << "==";
    
    Rcpp::List data = boot_data(y_0, n, K, t_size, beta_sub, cor, marginal, dispersion);
    const arma::vec &y = data["response"];
    const arma::mat &xx = data["covariates"];
    
    const arma::vec y_minus1 = y - 1;
    vector<double> lower_bd = bound(xx, y_minus1, beta_sub, marginal, dispersion), 
                   upper_bd = bound(xx, y, beta_sub, marginal, dispersion);
    
    vector<vector<double> >  upper_bdd = get_bd_deriv_arma_sub(xx, y, true, beta_sub, v_main, marginal, dispersion),
                             lower_bdd = get_bd_deriv_arma_sub(xx, y, false, beta_sub, v_main, marginal, dispersion);
    
    arma::vec tmp_score(p_sub, arma::fill::zeros);
    score_hessian(tmp_score, H, rho_v, v_rho, p_main_sub, p_sub, lower_bd, upper_bd, upper_bdd, lower_bdd, labeled_pairs, marginal);
    J += tmp_score * tmp_score.t();
  }
  H /= B;  J /= B; 
  
  if(Message_prog) Rcpp::Rcout << std::endl;
  
  arma::mat tmp_mat = arma::solve(H, J);
  
  arma::mat se_mat = tmp_mat * arma::inv(H);
  se = arma::diagvec(se_mat);
  
  return arma::trace(tmp_mat);
}
