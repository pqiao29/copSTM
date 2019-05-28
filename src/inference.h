#ifndef _INFERENCECPP_H
#define _INFERENCECPP_H
 
#include <RcppArmadillo.h>

double pois_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta
                        , arma::rowvec& score, arma::mat& hessian); 

double pois_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta, 
                 const arma::uvec& v, arma::rowvec& score, arma::mat& hessian);

double nbinom_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta,
                        arma::rowvec& score, arma::mat& hessian); 

double nbinom_infc(const arma::mat& xx, const arma::vec& y, const arma::vec& theta,
                   const arma::uvec& v, arma::rowvec& score, arma::mat& hessian);

#endif // _INFERENCECPP_H
