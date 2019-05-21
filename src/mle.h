
#ifndef _MLECPP_H
#define _MLECPP_H

#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <vector>
#include <map>

double mle(const arma::mat& xx, const arma::vec& y,
           arma::vec& beta, std::vector<double>& rho_v,
           int marginal, double& dispersion,
           const std::multimap<int, std::vector<int> >& labeled_pairs,
           int maxit, double eps);

double mle(const arma::mat& xx, const arma::vec& y,
           arma::vec& beta, std::vector<double>& rho_v,
           const std::multimap<int, std::vector<int> >& labeled_pairs,
           int maxit, double eps);


#endif // _MLECPP_H
