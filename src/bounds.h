
#ifndef _BOUNDSCPP_H
#define _BOUNDSCPP_H
 
#include <RcppArmadillo.h>
#include <vector>

std::vector<double> bound(const arma::mat& xx, const arma::vec& y,
                          const arma::vec& beta);


#endif // _BOUNDSCPP_H
