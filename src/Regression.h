#ifndef _REGRESSIONCPP_H
#define _REGRESSIONCPP_H
 
#include <RcppArmadillo.h>


Rcpp::List Regression(const arma::mat& xx, const arma::vec& y, int marginal, int maxit);

Rcpp::List Regression(const arma::mat& xx, const arma::vec& y, arma::vec theta, const arma::uvec& v_k, int marginal, int maxit);

#endif // _REGRESSIONCPP_H
