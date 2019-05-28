#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// Poisson, full model 
double pois_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta
                        , arma::rowvec& score, arma::mat& hessian) 
{
  arma::vec lmd = exp(xx * theta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    lik -= R::dpois(y(i), lmd(i), true);
    arma::rowvec tmp_row = xx.row(i);
    score -= (y(i) - lmd(i)) * tmp_row;
    hessian += lmd(i) * (tmp_row.t() * tmp_row);
  }
  return lik;
}

// Poisson, sub model
double pois_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta, 
                 const arma::uvec& v, arma::rowvec& score, arma::mat& hessian) 
{
  arma::vec lmd = exp(xx * theta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    lik -= R::dpois(y(i), lmd(i), true);
    arma::vec tmp_row = xx.row(i).t();  
              tmp_row = tmp_row.elem(find(v));
    score -= (y(i) - lmd(i)) * tmp_row.t();
    hessian += lmd(i) * (tmp_row * tmp_row.t());
  }
  return lik;
}

// Negative binomial, full model
double nbinom_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta,
                   arma::rowvec& score, arma::mat& hessian) 
{
  int p = theta.size();
  arma::vec beta = theta.head(p - 1);
  double overdispersion = theta(p - 1); 
  arma::vec mu = exp(xx * beta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    double y_i = y(i), mu_i = mu(i);
    lik -= R::dnbinom_mu(y_i, overdispersion, mu_i, true);
    arma::rowvec tmp_row = xx.row(i), tmp_score(p);
    
    // score, hessian
    tmp_score.head(p - 1) = tmp_row*overdispersion*(y_i - mu_i)/(overdispersion + mu_i);
    tmp_score(p - 1) = R::digamma(y_i + overdispersion) - 
                        R::digamma(overdispersion) + 
                        log(overdispersion/(overdispersion + mu_i)) + 
                        (mu_i - y_i)/(mu_i + overdispersion);

    score -= tmp_score;
    hessian += tmp_score.t() * tmp_score;
  }
  return lik;
}

// Negative binomial, sub model
double nbinom_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta,
                   const arma::uvec& v, arma::rowvec& score, arma::mat& hessian) 
{
  int p = theta.size(), p_v = sum(v);
  arma::uvec v_inuse = v.head(p - 1);
  arma::vec beta = theta.head(p - 1);
  double overdispersion = theta(p - 1); 
  
  arma::vec mu = exp(xx * beta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  
  for(int i = 0; i != y.size(); ++i){
    
    double y_i = y(i), mu_i = mu(i);
    
    lik -= R::dnbinom_mu(y_i, overdispersion, mu_i, true);
    
    arma::vec tmp_row = xx.row(i).t();  
              tmp_row = tmp_row.elem(find(v_inuse));
    
    arma::rowvec  tmp_score(p_v);
    // score, hessian
    tmp_score.head(p_v - 1) = tmp_row.t()*overdispersion*(y_i - mu_i)/(overdispersion + mu_i);
    tmp_score(p_v - 1) = R::digamma(y_i + overdispersion) - 
                         R::digamma(overdispersion) + 
                         log(overdispersion/(overdispersion + mu_i)) + 
                         (mu_i - y_i)/(mu_i + overdispersion);
    
    score -= tmp_score;
    hessian += tmp_score.t() * tmp_score;
  }
  
  return lik;
}

