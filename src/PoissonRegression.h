
#ifndef _POISSONREGRESSIONCPP_H
#define _POISSONREGRESSIONCPP_H
 
#include <RcppArmadillo.h>
 
void get_hessian(const arma::mat& xx, const arma::vec& y, const arma::vec& theta, arma::mat& hessian);
void get_hessian(const arma::mat& xx, const arma::vec& y, const arma::vec& theta, 
                 const arma::uvec& v, arma::mat& hessian);

double Poisson_Newton(const arma::mat& xx, const arma::vec& y, arma::vec& theta, bool happy = true);
 
double Poisson_Newton(const arma::mat& xx, const arma::vec& y, const arma::vec& theta, const arma::uvec& v, bool happy = true);
double Poisson_Newton(const arma::mat& xx, const arma::vec& y, arma::vec& theta, const arma::uvec& v, bool happy = true);


Rcpp::List PoissonRegression(const arma::mat& xx, const arma::vec& y);

#endif // _POISSONREGRESSIONCPP_H
