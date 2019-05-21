#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double bd_unary(const arma::vec& xx_row, const double y_i, bool upper,
                const arma::vec& beta, 
                int marginal, double dispersion){
  double lmd = exp(dot(xx_row, beta));
  
  if(upper){
    double tmp;
    if(marginal == 1){
      tmp = R::ppois(y_i, lmd, false, false);
    }else{
      tmp = R::pnbinom_mu(y_i, dispersion, lmd, false, false);
    }
    return -R::qnorm(tmp, 0, 1, true, false);
  }else{
    double tmp;
    if(marginal == 1){
      tmp = R::ppois((y_i - 1), lmd, false, false);
    }else{
      tmp = R::pnbinom_mu((y_i - 1), dispersion, lmd, false, false);
    }
    
    return -R::qnorm(tmp, 0, 1, true, false);
  }
}

#include <vector>
using std::vector;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> bd_deriv_unary_arma(const arma::vec& xx_row, const double y_i, 
                                   bool upper, const arma::vec& beta, 
                                   int marginal, double dispersion){
  
  if(std::isfinite(bd_unary(xx_row, y_i, upper, beta, marginal, dispersion))){
   
   std::vector<double> ret;
    arma::vec::size_type n = beta.size();
    double heps = 6.055454e-06;
    
    double tmp_ele;
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
      
      arma::vec tmp1 = beta, tmp2 = beta; 
      tmp1(ni) += heps;
      tmp2(ni) -= heps;
      
      tmp_ele = (bd_unary(xx_row, y_i, upper, tmp1, marginal, dispersion) - 
                 bd_unary(xx_row, y_i, upper, tmp2, marginal, dispersion))/(2 * heps);
      
      if(std::isfinite(tmp_ele))
        ret.push_back(tmp_ele);
      else 
        ret.push_back(0);
    }
    
    
    if(marginal == 2){
      double tmp1 = dispersion + heps, tmp2 = dispersion - heps;
      tmp_ele = (bd_unary(xx_row, y_i, upper, beta, marginal, tmp1) - 
                 bd_unary(xx_row, y_i, upper, beta, marginal, tmp2))/(2 * heps);
      if(std::isfinite(tmp_ele))
        ret.push_back(tmp_ele);
      else 
        ret.push_back(0);
    } 
    
    return ret;
  }else{
    std::vector<double> ret(beta.size(), 0);
    if(marginal == 2) ret.push_back(0);
    return ret;
  }
  
}


std::vector<std::vector<double> > get_bd_deriv_arma(const arma::mat& xx, const arma::vec& y, bool upper,
                                          const arma::vec& beta, 
                                          int marginal, double dispersion){
  std::vector<std::vector<double> > ret;
  arma::vec::size_type n = y.size();
  for(arma::vec::size_type ni = 0; ni != n; ++ni){
    arma::vec tmp_x = arma::conv_to<arma::vec>::from(xx.row(ni));
    ret.push_back(bd_deriv_unary_arma(tmp_x, y(ni), upper, beta, marginal, dispersion));
  }
  return ret;
}


/* ********************************************************************************
 * sub
 *
 ******************************************************************************** */

double bd_unary_arma_sub(const arma::vec& xx_row, const double y_i, bool upper,
                         const arma::vec& beta,
                         const arma::Col<int>& v_main, 
                         int marginal, double dispersion){
  
    double lmd = exp(dot(xx_row, (beta % v_main)));
    if(upper){
        double tmp; 
        if(marginal == 1){
          tmp = R::ppois(y_i, lmd, false, false);
        }else{
          tmp = R::pnbinom_mu(y_i, dispersion, lmd, false, false);
        } 
        return -R::qnorm(tmp, 0, 1, true, false);
    }else{
        double tmp; 
        if(marginal == 1){
          tmp = R::ppois((y_i - 1), lmd, false, false);
        }else{
          tmp = R::pnbinom_mu((y_i - 1), dispersion, lmd, false, false);
        }
      
        return -R::qnorm(tmp, 0, 1, true, false);
    }
}


vector<double> bd_deriv_unary_arma_sub(const arma::vec& xx_row, const double y_i, bool upper,
                                       const arma::vec& beta,
                                       const arma::Col<int>& v_main, 
                                       int marginal, double dispersion){
    
    if(std::isfinite(bd_unary_arma_sub(xx_row, y_i, upper, beta, v_main, marginal, dispersion))){
        vector<double> ret;
        arma::vec::size_type n = beta.size();
        double heps = 6.055454e-06;
        
    double tmp_ele;
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
            if(v_main(ni)){
                arma::vec tmp1 = beta, tmp2 = beta;
                tmp1(ni) += heps;
                tmp2(ni) -= heps;
                
                tmp_ele = (bd_unary_arma_sub(xx_row, y_i, upper, tmp1, v_main, marginal, dispersion) -
                           bd_unary_arma_sub(xx_row, y_i, upper, tmp2, v_main, marginal, dispersion))/(2 * heps);
                
                if(std::isfinite(tmp_ele))
                    ret.push_back(tmp_ele);
                else
                    ret.push_back(0);
            }
        }
    
    if(marginal == 2){
      double tmp1 = dispersion + heps, tmp2 = dispersion - heps;
      tmp_ele = (bd_unary_arma_sub(xx_row, y_i, upper, beta, v_main, marginal, tmp1) - 
                 bd_unary_arma_sub(xx_row, y_i, upper, beta, v_main, marginal, tmp2))/(2 * heps);
      if(std::isfinite(tmp_ele))
        ret.push_back(tmp_ele);
      else 
        ret.push_back(0);
    } 
    
    return ret;
    }else{
       int p = sum(v_main);
       if(marginal == 2) ++p; 
        vector<double> ret(p, 0);
        return ret;
    }
}

vector<vector<double> > get_bd_deriv_arma_sub(const arma::mat& xx, const arma::vec& y, bool upper,
                                              const arma::vec& beta,
                                              const arma::Col<int>& v_main, 
                                              int marginal = 1, double dispersion = 1){
    vector<vector<double> > ret;
    arma::vec::size_type n = y.size();
    for(arma::vec::size_type ni = 0; ni != n; ++ni){
        arma::vec tmp_x = arma::conv_to<arma::vec>::from(xx.row(ni));
        ret.push_back(bd_deriv_unary_arma_sub(tmp_x, y(ni), upper, beta, v_main, marginal, dispersion));
    }
    return ret;
}

