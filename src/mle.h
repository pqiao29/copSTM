
#ifndef _MLECPP_H
#define _MLECPP_H
 
#include <RcppArmadillo.h>
#include <vector>
#include <map>

double mle(const arma::mat& xx, const arma::vec& y,
           arma::vec& beta, std::vector<double>& rho_v,
           const std::multimap<int, std::vector<int> >& labeled_pairs,
           int maxit, double eps);

#endif // _MLECPP_H
