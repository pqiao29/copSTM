#include <Rcpp.h>
#include <vector>

#include "bvnorm.h"

double f_cond_cpp(const double& up, const double& down, const double& condi, const double& rho){
  if(std::isinf(condi)){
    return 0; 
  }else{
    return  R::dnorm(condi, 0, 1, false) 
    * (R::pnorm(up, rho*condi, sqrt(1 - rho*rho), true, false) 
         - R::pnorm(down, rho*condi, sqrt(1 - rho*rho), true, false));
  }
}

double h_cpp(const double& x, const double& y, const double& rho){
  if(std::isinf(x) || std::isinf(y)){
    return 0;
  }else{
    return exp(-(x*x - 2*rho*x*y + y*y)/(2*(1 - rho*rho)));
  }
}

/*
 *  Score for main effects
 */

double score_main_cpp(const std::vector<int>& pair_ind, const int& par_pos, const double& rho,
                      const std::vector<double>& lower_bd, const std::vector<double>& upper_bd, 
                      const std::vector<std::vector<double> >& upper_bdd, 
                      const std::vector<std::vector<double> >& lower_bdd){
  int i = pair_ind[0], j = pair_ind[1];
  double a1 = lower_bd[i], a2 = lower_bd[j], b1 = upper_bd[i], b2 = upper_bd[j];
  std::vector<double> bd1{a1, b1}, bd2{a2, b2};
  
  double comp_lik = dbnorm(bd1, bd2, rho);
  if(!comp_lik){
    return 0; 
  }else{
    double t1_up = f_cond_cpp(b2, a2, b1, rho),
      t1_down = f_cond_cpp(b2, a2, a1, rho),
      t2_up = f_cond_cpp(b1, a1, b2, rho),
      t2_down = f_cond_cpp(b1, a1, a2, rho);
    
    double b1_d = upper_bdd[i][par_pos], 
                              a1_d = lower_bdd[i][par_pos], 
                                                 b2_d = upper_bdd[j][par_pos], 
                                                                    a2_d = lower_bdd[j][par_pos];
    
    return  (b1_d*t1_up - a1_d*t1_down + b2_d*t2_up - a2_d*t2_down)/(comp_lik);
  }
}

/*
 *  Score for correlation parameters
 */
double score_rho_cpp(const std::vector<int>& pair_ind, double rho, 
                     const std::vector<double>& lower_bd, const std::vector<double>& upper_bd, 
                     const double& pi){
  int i = pair_ind[0], j = pair_ind[1];
  double a1 = lower_bd[i], a2 = lower_bd[j], b1 = upper_bd[i], b2 = upper_bd[j];
  std::vector<double> bd1{a1, b1}, bd2{a2, b2};
  
  double t1 = h_cpp(a1, b2, rho) + h_cpp(a2, b1, rho) - h_cpp(a1, a2, rho) - h_cpp(b1, b2, rho),
    t2 = dbnorm(bd1, bd2, rho);
  
  if(!t2){
    return 0; 
  }else{
    return -t1/(sqrt(1 - rho*rho)*2*pi*t2);
  }
}