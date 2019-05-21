
#ifndef _SCORE_HESSIANCPP_H
#define _SCORE_HESSIANCPP_H
 
#include <RcppArmadillo.h>
#include <map>
#include <vector>


void score_hessian(arma::vec& score, arma::mat& hessian,
                   const std::vector<double>& rho_v,
                   const std::vector<double>& lower_bd, const std::vector<double>& upper_bd,
                   const std::vector<std::vector<double> >& upper_bdd,
                   const std::vector<std::vector<double> >& lower_bdd,
                   const std::multimap<int, std::vector<int> >& labeled_pairs, 
                   int marginal = 1);

void score_hessian(arma::vec& score, arma::mat& hessian,
                   const std::vector<double>& rho_v,  const arma::Col<int>& v_rho,
                   const int& p_main, const int& p,
                   const std::vector<double>& lower_bd, const std::vector<double>& upper_bd,
                   const std::vector<std::vector<double> >& upper_bdd,
                   const std::vector<std::vector<double> >& lower_bdd,
                   const std::multimap<int, std::vector<int> >& labeled_pairs, 
                   int marginal = 1);

double get_std_err(const arma::mat& xx, const arma::vec& y, 
                   const arma::vec& beta, const std::vector<double>& rho_v, 
                   int t_size, int d, arma::vec& se, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs0, 
                   int marginal = 1, double dispersion = 1);

double get_std_err(const arma::mat& xx, const arma::vec& y, 
                   const arma::vec& beta, std::vector<double> rho_v, const int t_size, 
                   const arma::Col<int>& v_main, const arma::Col<int>& v_rho, 
                   int d, arma::vec& se, 
                   const std::multimap<int, std::vector<int> >& labeled_pairs0, 
                   int marginal = 1, double dispersion = 1);

#endif // _SCORE_HESSIANCPP_H
