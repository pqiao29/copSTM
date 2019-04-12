#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <iterator>
using namespace std; 

#include <algorithm>
using std::transform;

#include <Rcpp.h>
using namespace Rcpp;

inline double f_exp(double x){ return exp(x); }
inline double f_ppois(double x, double lambda){ return R::ppois(x, lambda, false, false); }
inline double f_qnorm(double x){return -R::qnorm(x, 0, 1, true, false);}


vector<double> bound(const arma::mat& xx, const arma::vec& y, const arma::vec& beta){
  vector<double> ret = arma::conv_to< vector<double> >::from( xx * beta );
  transform(ret.begin(), ret.end(), ret.begin(), f_exp);
  transform(y.cbegin(), y.cend(), ret.begin(), ret.begin(), f_ppois);
  transform(ret.begin(), ret.end(), ret.begin(), f_qnorm);
  return ret;
}

