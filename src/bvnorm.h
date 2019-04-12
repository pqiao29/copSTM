
#ifndef _BVNORMCPP_H
#define _BVNORMCPP_H

#include <vector>

double Phi(double value);

double dbnorm_inf( const double& x,
             const double& y,
             const double& rho );

double dbnorm( const std::vector<double>& x,
              const std::vector<double>& y,
               const double& rho );

#endif // _BVNORMCPP_H
