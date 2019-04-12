
#ifndef _NEAR_PDCPP_H
#define _NEAR_PDCPP_H
 
#include <RcppArmadillo.h>

void nearPD_cpp(arma::mat& X, bool corr = false
                , double eig_tol   = 1e-6
                , double conv_tol  = 1e-7
                , double posd_tol  = 1e-8
                , int maxit = 200
                );

#endif // _NEAR_PDCPP_H
