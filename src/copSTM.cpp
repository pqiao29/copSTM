#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Regression.h"
#include "inference.h"
#include "bounds.h"
#include "bounds_deriv.h"
#include "labeled_pairs.h"
#include "score_hessian.h"
#include "lik.h"
#include "mle.h"
#include "boot_penalty.h"
#include "make_cor.h"

#include <vector>
#include <map>

// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
Rcpp::List copSTM_cpp(const arma::mat& x, const arma::vec& y, 
                      int marginal,
                      const bool temporal, const int cor_type,
                      const int K, const int n, int maxit, double eps, 
                      bool std_err, int B = 0, bool Message_prog = false){
  
  // Input xx need to include the 1 column for intercept
  
  // Initialize
  auto ini = Regression(x, y, marginal, maxit);
  arma::vec beta_ini = ini["paramters"];
  double lik = ini["likelihood"];
  if(std::isnan(lik) || std::isinf(lik)) throw Rcpp::exception("Unsuccessful glm fit.", false);
  
  if(cor_type == 4){
    
    int p = x.n_cols;
    if(marginal == 2) ++p;
    
    arma::mat hessian(p, p); arma::rowvec score(p);
    
    if(marginal == 1) pois_infc(x, y, beta_ini, score, hessian);
    if(marginal == 2) nbinom_infc(x, y, beta_ini, score, hessian);
    
    arma::vec se = arma::diagvec(arma::inv_sympd(hessian));
    
    return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                              Rcpp::Named("main") = beta_ini.head(p - 1), 
                              Rcpp::Named("dispersion") = beta_ini(p - 1),
                              Rcpp::Named("std_err") = se);
  }
  
  //labeled_pairs
  int d = K*n*n; 
  int t_size = y.size()/d, p_rho;
  const std::multimap<int, std::vector<int> > labeled_pairs = get_pairs(K, n, p_rho, t_size, cor_type); //override p_rho
  std::vector<double> rho_v(p_rho, 0.0);
  
  // Copula estimation
  arma::vec beta;
  double dispersion = 1;
  if( marginal == 1){
    beta = beta_ini;
    //lik = mle(x, y, beta_ini, rho_v, labeled_pairs, maxit, eps);
  }else{
    beta = beta_ini.head(x.n_cols);
    dispersion = arma::conv_to<double>::from(beta_ini.tail(1));
  }
  lik = mle(x, y, beta, rho_v, marginal, dispersion, labeled_pairs, maxit, eps);
   // override only beta and rho_v
 
  
  if(std_err){ 
    
    arma::vec se;
    int p_main = x.n_cols, p = p_main + p_rho;
    const std::multimap<int, std::vector<int> > labeled_pairs0 = get_pairs(K, n, p_rho, 1, cor_type); 
    
    if(temporal){ // Standard error for temporal setting
      arma::vec y_0 = y.head(d);
      arma::mat corr = cor_mat(rho_v, labeled_pairs0, d);
      if(marginal == 2) ++p; 
      boot_CLIC_penalty(y_0, n, K, t_size, beta, se, corr, 
                        rho_v, p_main, p, marginal, dispersion, labeled_pairs, B, Message_prog);
    }else{
      get_std_err(x, y, beta, rho_v, t_size, d, se,
                  labeled_pairs0, marginal, dispersion);
    }
    
    return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                              Rcpp::Named("main") = beta, 
                              Rcpp::Named("rho") = rho_v, 
                              Rcpp::Named("dispersion") = dispersion,
                              Rcpp::Named("std_err") = se);
 
  }
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("main") = beta, 
                            Rcpp::Named("rho") = rho_v, 
                            Rcpp::Named("dispersion") = dispersion);
  
}