
#ifndef _MLE_SUBCPP_H
#define _MLE_SUBCPP_H
 
#include <RcppArmadillo.h>
#include <vector>
#include <map>


Rcpp::List mle_sub(double& l, const arma::mat& xx, const arma::vec& y,
                  arma::vec beta, std::vector<double> rho_v,
                  const arma::Col<int>& v_main, const arma::Col<int>& v_rho,
                  const std::multimap<int, std::vector<int> >& labeled_pairs,
                  int maxit, double eps);


#endif // _MLE_SUBCPP_H
