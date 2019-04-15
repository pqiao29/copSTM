#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "grid_ring.h"
#include "PoissonRegression.h"

// [[Rcpp::plugins(cpp11)]]


Rcpp::List data_indpt(const arma::mat& dat, const int n){
  
  std::map<int, std::vector<int> > nbr = grid_ring(n, 1);
  
  int t_size = dat.col(0).max() + 1, K = dat.col(1).max(), n_lattice = n*n;
  
  // count
  arma::cube obs_cnt(K, n_lattice, t_size, arma::fill::zeros);
  for(int i = 0; i != dat.n_rows; ++i){
    ++obs_cnt((dat(i, 1) - 1), (dat(i, 2) - 1), dat(i, 0));
  }
  
  arma::mat response((n_lattice * (t_size - 1)), K); //uninitialized
  arma::mat covariate((n_lattice * (t_size - 1)), K + 1); covariate.col(0).ones();
  int ind_cov = 0, ind_res = 0; 
  for(int t = 1; t != t_size; ++t){
    // response
    response.rows(ind_res, (ind_res + n_lattice - 1)) = obs_cnt.slice(t).t();
    ind_res += n_lattice;
    // covariate
    for(int i = 0; i != n_lattice; ++i){
      std::vector<int> tmp_nb = nbr[i + 1];
      double tmp; 
      for(int k = 0; k != K; ++k){
        tmp = 0; 
        std::vector<int>::const_iterator iter = tmp_nb.cbegin();
        while(iter != tmp_nb.cend()){
          tmp += log(obs_cnt(k, (*iter++ - 1), t - 1) + 1);
        }
        tmp /= tmp_nb.size();
        covariate(ind_cov, k + 1) = tmp; 
      }
      ++ind_cov;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("covariates") = covariate, 
                            Rcpp::Named("T") = t_size, 
                            Rcpp::Named("K") = K);
  
} 

// [[Rcpp::export]]
Rcpp::List idptSTM_cpp(const arma::mat& dat, const int n_lattice, const int maxit){
  
  /*
   * Requirement of data columns: Timepoint (0, 1, 2, ...), group (1, 2, 3, ...), tile (1, 2, ... n*n)
   */
  
  auto in_data = data_indpt(dat, n_lattice);
  arma::mat x = in_data["covariates"]; 
  arma::mat y = in_data["response"];
  int K = in_data["K"];
  double lik = 0; 
  
  arma::rowvec beta0(K); arma::mat beta(K, K); arma::mat se((K + 1), K);
  int ind_b0 = 0;
  for(int k = 0; k != K; ++k){
    int tmp_maxit = maxit; bool happy = true;
    // est
    int p = x.n_cols;
    arma::vec theta(p); theta.zeros();
    theta(0) = sum(y.col(k))/y.size();
    
    double l = Poisson_Newton(x, y.col(k), theta, tmp_maxit, happy);
    //if(!happy) throw Rcpp::exception("Unsuccessful glm fit.", false);

    beta0(ind_b0++) = theta(0);
    beta.col(k) = theta.tail(K);
    lik += l;
    
    // std_err 
    arma::mat hessian(p, p);
    get_hessian(x, y.col(k), theta, hessian);
    se.col(k) = arma::diagvec(arma::inv_sympd(hessian));
  }

  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("intercept") = beta0, 
                            Rcpp::Named("main_effects") = beta, 
                            Rcpp::Named("se") = se);
  
}
  