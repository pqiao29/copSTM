
#ifndef _BOOT_DATACPP_H
#define _BOOT_DATACPP_H
 
#include <RcppArmadillo.h>

Rcpp::List boot_data(const arma::vec& y_0,
                     const int n, const int K, const int t_size,
                     const arma::vec& beta, const arma::mat& cor);

Rcpp::List gen_data( int t_size, int d, const arma::vec& beta,
                     const arma::mat& cor);


#endif // _BOOT_DATACPP_H
