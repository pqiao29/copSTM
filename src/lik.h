
#ifndef _LIKCPP_H
#define _LIKCPP_H
 
#include <Rcpp.h>
#include <vector>
#include <map>
using namespace Rcpp;

double lik(const std::vector<double>& rho_v,
           const std::multimap<int, std::vector<int> >& labeled_pairs,
           const std::vector<double>& lower_bd,
           const std::vector<double>& upper_bd);

#endif // _LIKCPP_H
