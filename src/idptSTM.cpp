#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "grid_ring.h"
#include "Regression.h"
#include "inference.h"

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
Rcpp::List idptSTM_cpp(const arma::mat& dat, const int n, const int marginal, 
                       const int maxit, bool fit_plot){
  
  /*
   * Requirement of data columns: Timepoint (0, 1, 2, ...), group (1, 2, 3, ...), tile (1, 2, ... n*n)
   */
  
  auto in_data = data_indpt(dat, n);
  arma::mat x = in_data["covariates"]; 
  arma::mat y = in_data["response"];
  
  int K = in_data["K"];
  double lik = 0; 
  int t_size = in_data["T"];
  arma::Mat<int> fitted(K, t_size - 1), obsved(K, t_size - 1);
  
  arma::rowvec beta0(K), dispersion(K); 
  arma::mat beta(K, K), se((K + 1), K);
  arma::rowvec se_dispersion(K);
  
  int ind_b0 = 0;
  for(int k = 0; k != K; ++k){
    int tmp_maxit = maxit;
     
    // est
    auto ini = Regression(x, y.col(k), marginal, tmp_maxit);
    arma::vec theta = ini["paramters"];
    double l = ini["likelihood"];
    lik += l;
    //double l = Poisson_Newton(x, y.col(k), theta, tmp_maxit, happy);
    //if(!happy) throw Rcpp::exception("Unsuccessful glm fit.", false);

    beta0(ind_b0++) = theta(0);
    beta.col(k) = theta.subvec(1, K);
    if(marginal == 2) dispersion(k) = theta(K + 1);
    
    // std_err 
    int p; 
    if(marginal == 1)  p = x.n_cols;
    if(marginal == 2)  p = x.n_cols + 1;
    
    arma::mat hessian(p, p);
    arma::rowvec score(p);
    if(marginal == 1) pois_infc(x, y.col(k), theta, score, hessian);
    if(marginal == 2) nbinom_infc(x, y.col(k), theta, score, hessian);
    
    //get_hessian(x, y.col(k), theta, hessian);
    if(marginal == 1) se.col(k) = arma::diagvec(arma::inv_sympd(hessian));
    if(marginal == 2){
      arma::vec tmp_se = arma::diagvec(arma::inv_sympd(hessian));
      se.col(k) = tmp_se.head(K + 1);
      se_dispersion(k) = arma::conv_to<double>::from(tmp_se.tail(1));
    } 
    
    // for goodness-of-fit curves
    if(fit_plot && marginal == 1){
      int n_lattice = n*n; 
      arma::vec fit_lmd = exp(x * theta);
      fit_lmd.for_each( [](arma::vec::elem_type& l){ l = R::rpois(l); } );
      for(int t = 0; t != t_size - 1; ++t){
        fitted(k, t) = sum(fit_lmd.subvec(t*n_lattice, ((t + 1)*n_lattice - 1)));
        obsved(k, t) = sum((y.col(k)).subvec(t*n_lattice, ((t + 1)*n_lattice - 1)));
      }
    }
  }

  if(fit_plot && marginal == 1){
    return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                              Rcpp::Named("intercept") = beta0, 
                              Rcpp::Named("main_effects") = beta, 
                              Rcpp::Named("se") = se, 
                              Rcpp::Named("fit") = fitted, 
                              Rcpp::Named("obs") = obsved);
  }else{
    
    if(marginal == 2)
      return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                                Rcpp::Named("intercept") = beta0, 
                                Rcpp::Named("main_effects") = beta, 
                                Rcpp::Named("dispersion") = dispersion, 
                                Rcpp::Named("se") = se, 
                                Rcpp::Named("se_dispersion") = se_dispersion);
      
  }
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("intercept") = beta0, 
                            Rcpp::Named("main_effects") = beta, 
                            Rcpp::Named("se") = se);
  
}
  