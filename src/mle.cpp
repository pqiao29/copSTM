#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "bounds.h"
#include "bounds_deriv.h"
#include "score_hessian.h"
#include "lik.h"


double mle(const arma::mat& xx, const arma::vec& y, 
               arma::vec& beta, std::vector<double>& rho_v, 
               int marginal, double& dispersion,
               const std::multimap<int, std::vector<int> >& labeled_pairs,
               int maxit, double eps){
  
  int p_main = beta.size(), p_rho = rho_v.size(), p = p_main + p_rho;
  if(marginal == 2) ++p;
  
  // theta (concatenaing beta, rho_v, dispersion)
  arma::vec theta(p);
  theta.head(p_main) = beta;
  theta.subvec(p_main, (p_main + p_rho - 1)) = arma::vec(rho_v);
  if(marginal == 2) theta(p - 1) = dispersion;
  
  bool converged = false; 
  double l_prev = 0; 
  arma::vec score(p); arma::mat hessian(p, p); 
  int mle_iteration = 0; 
  
  // Workhorse
  while(!converged && mle_iteration++ != maxit){
    
    // Bounds for integrals in Gaussian cdf
    const arma::vec y_minus1 = y - 1;
    std::vector<double> lower_bd = bound(xx, y_minus1, beta, marginal, dispersion), 
                        upper_bd = bound(xx, y, beta, marginal, dispersion);
    
    std::vector<std::vector<double> > upper_bdd = get_bd_deriv_arma(xx, y, true, beta, marginal, dispersion),
                                      lower_bdd = get_bd_deriv_arma(xx, y, false, beta, marginal, dispersion);
    
    // Update
    score.zeros(); hessian.zeros();
    score_hessian(score, hessian, rho_v, lower_bd, upper_bd, 
                  upper_bdd, lower_bdd, labeled_pairs, marginal);
    
    if(marginal == 1){
      theta += arma::solve(hessian, score);
    }else{
      theta += arma::inv(hessian) * arma::vec(score);
    }
    
    beta = theta.head(p_main);
    rho_v = arma::conv_to< std::vector<double> >::from(theta.subvec(p_main, (p_main + p_rho - 1)));  
    if(marginal == 2) dispersion = theta(p - 1);
    
    if(dispersion <= 0) dispersion= 1;

    // Check convergence
    double l = lik(rho_v, labeled_pairs, lower_bd, upper_bd);
    
    if(!std::isnan(l) && !std::isinf(l)){
      if(std::abs(l - l_prev) < eps){
        converged = true;
      }
      l_prev = l; 
    }else{
      throw Rcpp::exception("undefined likelihood produced", false);
    }
    
  }
  
  if(!converged) Rcpp::warning("Maximum likelihood estimation not converged");
  
  return l_prev;
}

