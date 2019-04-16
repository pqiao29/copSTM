#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bounds.h"
#include "bounds_deriv.h"
#include "score_hessian.h"
#include "lik.h"

#include <vector>
#include <map>

// [[Rcpp::plugins(cpp11)]]
Rcpp::List mle_sub(double& l, const arma::mat& xx, const arma::vec& y, 
                   arma::vec beta, std::vector<double> rho_v, 
                   const arma::Col<int>& v_main, const arma::Col<int>& v_rho,
                   const std::multimap<int, std::vector<int> >& labeled_pairs,
                   int maxit, double eps){ // override l
  
  int p_main = sum(v_main), p = p_main + sum(v_rho);
  arma::vec beta_sub = beta % v_main;
  // theta (concatenaing beta0, beta, rho_v)
  arma::vec theta(nonzeros(beta_sub));
  arma::vec rho_zeros = arma::vec(rho_v) % v_rho;
  rho_v = arma::conv_to<std::vector<double> >::from(rho_zeros);
  theta = arma::join_cols(theta, nonzeros(rho_zeros));

  bool converged = false; 
  double l_prev = 0; 
  arma::vec score(p);  arma::mat hessian(p, p); 
  int iteration = 0; 
  
  while(!converged && iteration++ != maxit){
    // Bounds for integrals in Gaussian cdf
    const arma::vec y_minus1 = y - 1;
    std::vector<double> lower_bd = bound(xx, y_minus1, beta_sub), 
                        upper_bd = bound(xx, y, beta_sub);
    
    std::vector<std::vector<double> > upper_bdd = get_bd_deriv_arma_sub(xx, y, true, beta_sub, v_main),
                                      lower_bdd = get_bd_deriv_arma_sub(xx, y, false, beta_sub, v_main);
    
    // Update
    score.zeros(); hessian.zeros();
    score_hessian(score, hessian, rho_v, v_rho, p_main, p, lower_bd, upper_bd, upper_bdd, lower_bdd, labeled_pairs);
    
    //theta += arma::inv(hessian) * arma::vec(score);
    theta += arma::solve(hessian, score);
    // theta -> beta_sub, rho_v
    arma::vec::const_iterator theta_ptr = theta.cbegin();
    arma::Col<int>::const_iterator main_ptr = v_main.cbegin(), 
                                   rho_ptr = v_rho.cbegin();
    arma::vec::iterator beta_ptr = beta_sub.begin(); 
    std::vector<double>::iterator rho_v_ptr = rho_v.begin();
    
    while(main_ptr != v_main.cend()){
      if(*main_ptr++)
       *beta_ptr = *theta_ptr++;
      ++beta_ptr;
    }
    while(rho_ptr != v_rho.cend()){
      if(*rho_ptr++)
        *rho_v_ptr = *theta_ptr++;
      ++rho_v_ptr;
    }
     
    // Check convergence
    l = lik(rho_v, labeled_pairs, lower_bd, upper_bd);
    if(!std::isnan(l) && !std::isinf(l)){
      if(std::abs(l - l_prev) < eps){
        converged = true;
      }
      l_prev = l; 
    }else break;
  
  }
  
  // No warning generated if not converged. 
  // If full model estimation is successful while a submodel is not, then this submodel should not be selected
  
  return Rcpp::List::create(Rcpp::Named("beta") = beta_sub,
                            Rcpp::Named("rho") = rho_v );
}

