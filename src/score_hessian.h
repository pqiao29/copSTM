
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
                   const std::multimap<int, std::vector<int> >& labeled_pairs);

void score_hessian(arma::vec& score, arma::mat& hessian,
                   const std::vector<double>& rho_v,  const arma::Col<int>& v_rho,
                   const int& p_main, const int& p,
                   const std::vector<double>& lower_bd, const std::vector<double>& upper_bd,
                   const std::vector<std::vector<double> >& upper_bdd,
                   const std::vector<std::vector<double> >& lower_bdd,
                   const std::multimap<int, std::vector<int> >& labeled_pairs);

#endif // _SCORE_HESSIANCPP_H