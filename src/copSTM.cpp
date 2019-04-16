#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "PoissonRegression.h"
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
Rcpp::List copSTM_cpp(const arma::mat& x, const arma::vec& y, const int cor_type,
                      int K, int n, int maxit, double eps, 
                      bool std_err, int B = 0, bool Message_prog = false){
  
  // Input xx need to include the 1 column for intercept
  
  int d = K*n*n; 
  // Initialize
  int tmp_maxit = maxit;
  auto ini = PoissonRegression(x, y, maxit);
  arma::vec beta_ini = ini["paramters"];
  if(!ini["happyending"]) throw Rcpp::exception("Unsuccessful glm fit.", false);
  if(Message_prog) Rcpp::Rcout << "glm initial converged in " << tmp_maxit - maxit << " iterations" << std::endl;
  
  //labeled_pairs
  int t_size = y.size()/d, p_rho;
  const std::multimap<int, std::vector<int> > labeled_pairs = get_pairs(K, n, p_rho, t_size, cor_type); //override p_rho
  std::vector<double> rho_v(p_rho, 0.0);
    
  // Fit model
  double lik = mle(x, y, beta_ini, rho_v, labeled_pairs, maxit, eps); // override beta, rho_v
 
  // Standard error
  arma::vec se, y_0 = y.head(d);
  if(std_err){
    int p_main = x.n_cols, p = p_main + p_rho;
    const std::multimap<int, std::vector<int> > labeled_pairs0 = get_pairs(K, n, p_rho, 1, cor_type); 
    arma::mat corr = cor_mat(rho_v, labeled_pairs0, d);
    
    boot_CLIC_penalty(y_0, n, K, t_size, beta_ini, se, corr,   
                      rho_v, p_main, p, labeled_pairs, B, Message_prog);
    
    return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                              Rcpp::Named("main") = beta_ini, 
                              Rcpp::Named("rho") = rho_v, 
                              Rcpp::Named("std_err") = se);
  }
   
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("main") = beta_ini, 
                            Rcpp::Named("rho") = rho_v);

}