#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bounds.h"
#include "bounds_deriv.h"
#include "score_hessian.h"
#include "lik.h"

#include <vector>
#include <map>

// [[Rcpp::plugins(cpp11)]]
double mle(const arma::mat& xx, const arma::vec& y, 
               arma::vec& beta, std::vector<double>& rho_v, 
               const std::multimap<int, std::vector<int> >& labeled_pairs,
               int maxit, double eps){
  
  int p_main = beta.size(), p_rho = rho_v.size(), p = p_main + p_rho;
  // theta (concatenaing beta0, beta, rho_v)
  arma::vec theta(beta);
  theta = arma::join_cols(theta, arma::vec(rho_v));
  
  bool converged = false; 
  double l_prev = 0; 
  arma::vec score(p); arma::mat hessian(p, p); 
  int mle_iteration = 0; 
  while(!converged && mle_iteration++ != maxit){
    
    // Bounds for integrals in Gaussian cdf
    const arma::vec y_minus1 = y - 1;
    std::vector<double> lower_bd = bound(xx, y_minus1, beta), 
                        upper_bd = bound(xx, y, beta);
    
    std::vector<std::vector<double> > upper_bdd = get_bd_deriv_arma(xx, y, true, beta),
                                      lower_bdd = get_bd_deriv_arma(xx, y, false, beta);
    
    // Update
    score.zeros(); hessian.zeros();
    score_hessian(score, hessian, rho_v, lower_bd, upper_bd, upper_bdd, lower_bdd, labeled_pairs);
    
    //theta += arma::inv(hessian) * arma::vec(score);
    theta += arma::solve(hessian, score);
    beta = theta.head(p_main);
    rho_v = arma::conv_to< std::vector<double> >::from(theta.subvec(p_main, p - 1));  

    // Check convergence
    double l = lik(rho_v, labeled_pairs, lower_bd, upper_bd);
    if(!std::isnan(l) && !std::isinf(l)){
      if(std::abs(l - l_prev) < eps){
        converged = true;
      }
      l_prev = l; 
    }
  }
  
  if(!converged) Rcpp::warning("Maximum likelihood estimation not converged");
  
  return l_prev;
}
