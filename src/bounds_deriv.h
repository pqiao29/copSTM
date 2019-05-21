
#ifndef _BOUNDS_ARMACPP_H
#define _BOUNDS_ARMACPP_H
 
#include <RcppArmadillo.h>
#include <vector>

std::vector<std::vector<double> > get_bd_deriv_arma(const arma::mat& xx,
                                                    const arma::vec& y, bool upper,
                                                    const arma::vec& beta,
                                                    int marginal = 1, double dispersion = 1);

std::vector<std::vector<double> > get_bd_deriv_arma_sub(const arma::mat& xx,
                                                        const arma::vec& y, bool upper,
                                                        const arma::vec& beta,
                                                        const arma::Col<int>& v_main, 
                                                        int marginal = 1, double dispersion = 1);

#endif // _BOUNDS_ARMACPP_H
