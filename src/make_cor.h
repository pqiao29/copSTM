
#ifndef _MAKE_CORCPP_H
#define _MAKE_CORCPP_H
 
#include <RcppArmadillo.h>
#include <vector>
#include <map>


arma::mat cor_mat(const std::vector<double>& rho_v,
                  const std::multimap<int, std::vector<int> >& labeled_pairs, int d);

#endif // _MAKE_CORCPP_H
