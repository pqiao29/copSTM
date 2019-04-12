
#ifndef _BOOT_PENALTYCPP_H
#define _BOOT_PENALTYCPP_H
 
#include <RcppArmadillo.h>
#include <map>
#include <vector>

double boot_CLIC_penalty(  const arma::vec& y_0,
                           const int n, const int K, const int t_size, 
                           const arma::vec& beta,
                           const arma::mat& cor,   
                           const std::vector<double>& rho_v,
                           const int& p_main, const int& p,
                           const std::multimap<int, std::vector<int> >& labeled_pairs, 
                           int B, bool Message_prog);

double boot_CLIC_penalty(  const arma::vec& y_0,
                           const int n, const int K, const int t_size, 
                           const arma::vec& beta, arma::vec& se,
                           const arma::mat& cor,   
                           const std::vector<double>& rho_v,
                           const int& p_main, const int& p,
                           const std::multimap<int, std::vector<int> >& labeled_pairs, 
                           int B, bool Message_prog);


double boot_CLIC_penalty_sub(  const arma::vec& y_0,
                               const int n, const int K, const int t_size, 
                               const arma::vec& beta, 
                               const arma::mat& cor, 
                               const std::vector<double>& rho_v,
                               const arma::Col<int>& v_main, const arma::Col<int>& v_rho,
                               const int& p_main_sub, const int& p_sub,
                               const std::multimap<int, std::vector<int> >& labeled_pairs, 
                               int B, bool Message_prog);

double boot_CLIC_penalty_sub(  const arma::vec& y_0,
                           const int n, const int K, const int t_size, 
                           const arma::vec& beta, arma::vec& se,
                           const arma::mat& cor, 
                           const std::vector<double>& rho_v,
                           const arma::Col<int>& v_main, const arma::Col<int>& v_rho,
                           const int& p_main_sub, const int& p_sub,
                           const std::multimap<int, std::vector<int> >& labeled_pairs, 
                           int B, bool Message_prog);

#endif // _BOOT_PENALTYCPP_H
