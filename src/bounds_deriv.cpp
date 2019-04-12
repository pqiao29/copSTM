#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


double bd_unary_arma(const arma::vec& xx_row, const double y_i, bool upper,
                     const arma::vec& beta){
  double lmd = exp(dot(xx_row, beta));
  if(upper){
    double tmp = R::ppois(y_i, lmd, false, false);
    return -R::qnorm(tmp, 0, 1, true, false);
  }else{
    double tmp = R::ppois((y_i - 1), lmd, false, false);
    return -R::qnorm(tmp, 0, 1, true, false);
  }
}

#include <vector>
using std::vector;
// [[Rcpp::plugins(cpp11)]]
vector<double> bd_deriv_unary_arma(const arma::vec& xx_row, const double y_i, bool upper,
                              const arma::vec& beta){
  
  if(std::isfinite(bd_unary_arma(xx_row, y_i, upper, beta))){
    vector<double> ret;
    arma::vec::size_type n = beta.size();
    double heps = 6.055454e-06;
    
    double tmp_ele;
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
      
      arma::vec tmp1 = beta, tmp2 = beta; 
      tmp1(ni) += heps;
      tmp2(ni) -= heps;
      
      tmp_ele = (bd_unary_arma(xx_row, y_i, upper, tmp1) - bd_unary_arma(xx_row, y_i, upper, tmp2))/(2 * heps);
      
      if(std::isfinite(tmp_ele))
        ret.push_back(tmp_ele);
      else 
        ret.push_back(0);
    }
    return ret;
  }else{
    vector<double> ret(beta.size(), 0);
    return ret;
  }
}


vector<vector<double> > get_bd_deriv_arma(const arma::mat& xx, const arma::vec& y, bool upper,
                                          const arma::vec& beta){
  vector<vector<double> > ret;
  arma::vec::size_type n = y.size();
  for(arma::vec::size_type ni = 0; ni != n; ++ni){
    arma::vec tmp_x = arma::conv_to<arma::vec>::from(xx.row(ni));
    ret.push_back(bd_deriv_unary_arma(tmp_x, y(ni), upper, beta));
  }
  return ret;
}


/*
 * sub
 *
 */

double bd_unary_arma_sub(const arma::vec& xx_row, const double y_i, bool upper,
                         const arma::vec& beta,
                         const arma::Col<int>& v_main){
    double lmd = exp(dot(xx_row, (beta % v_main)));
    if(upper){
        double tmp = R::ppois(y_i, lmd, false, false);
        return -R::qnorm(tmp, 0, 1, true, false);
    }else{
        double tmp = R::ppois((y_i - 1), lmd, false, false);
        return -R::qnorm(tmp, 0, 1, true, false);
    }
}


vector<double> bd_deriv_unary_arma_sub(const arma::vec& xx_row, const double y_i, bool upper,
                                       const arma::vec& beta,
                                       const arma::Col<int>& v_main){
    
    if(std::isfinite(bd_unary_arma_sub(xx_row, y_i, upper, beta, v_main))){
        vector<double> ret;
        arma::vec::size_type n = beta.size();
        double heps = 6.055454e-06;
        
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
            if(v_main(ni)){
                arma::vec tmp1 = beta, tmp2 = beta;
                tmp1(ni) += heps;
                tmp2(ni) -= heps;
                
                double tmp_ele = (bd_unary_arma(xx_row, y_i, upper, tmp1) -
                           bd_unary_arma(xx_row, y_i, upper, tmp2))/(2 * heps);
                
                if(std::isfinite(tmp_ele))
                    ret.push_back(tmp_ele);
                else
                    ret.push_back(0);
            }
        }
        return ret;
    }else{
        vector<double> ret((sum(v_main)), 0);
        return ret;
    }
}

vector<vector<double> > get_bd_deriv_arma_sub(const arma::mat& xx, const arma::vec& y, bool upper,
                                              const arma::vec& beta,
                                              const arma::Col<int>& v_main){
    vector<vector<double> > ret;
    arma::vec::size_type n = y.size();
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
        arma::vec tmp_x = arma::conv_to<arma::vec>::from(xx.row(ni));
        ret.push_back(bd_deriv_unary_arma_sub(tmp_x, y(ni), upper, beta, v_main));
    }
    return ret;
}

