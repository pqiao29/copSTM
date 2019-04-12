
#ifndef _UTILITYCPP_H
#define _UTILITYCPP_H
 
#include <cmath>
#include <vector>

double f_cond_cpp(const double& up, const double& down, const double& condi, const double& rho);

double h_cpp(const double& x, const double& y, const double& rho);

double score_main_cpp(const std::vector<int>& pair_ind, const int& par_pos, const double& rho,
                      const std::vector<double>& lower_bd, const std::vector<double>& upper_bd,
                      const std::vector<std::vector<double> >& upper_bdd,
                      const std::vector<std::vector<double> >& lower_bdd);

double score_rho_cpp(const std::vector<int>& pair_ind, double rho,
                     const std::vector<double>& lower_bd, const std::vector<double>& upper_bd, const double& pi);

#endif // _UTILITYCPP_H
