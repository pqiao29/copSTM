#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 

#include"inference.h"


inline bool keep_going(double lcrit, double leps, int iteration,
                       int max_iterations, bool step_halving = false) {
  if (step_halving) {
    // If the last iteration used step halving we're not done.
    return true;
  } else if (iteration >= max_iterations) {
    // If the maximum number of iterations is exceeded it is time
    // to bail out.
    return false;
  } else if (lcrit > leps) {
    // If the convergence criterion exceeds epsilon then there is
    // more work to do.
    return true;
  } else {
    // Mission accomplished!
    return false;
  }
}

inline bool BAD(double lcrit, double epsilon) {
  return ((std::fabs(lcrit) > epsilon) && (lcrit < 0));
}


// full model 
double Newton(const mat& xx, const vec& y, vec& theta, 
              double inference(const mat&, const vec&, const vec&, rowvec&, mat&),
              int& max_iterations, bool& happy, int marginal = 1){
  //Initialize paramters 
  int p = theta.size();
  rowvec score(p);
  mat hessian(p, p);
  
  // Prepre for mle
  int iteration = 0;
  double lcrit = 1, lik = 0;
  double oldlik = inference(xx, y, theta, score, hessian);
  
  int step_halving = 0, total_step_halving = 0;
  const int max_step_halving = 10, max_total_step_halving = 50;
  
  // Newton update
  while(keep_going(lcrit, 1e-5, iteration, max_iterations)){
    
    if ( score.has_inf() || hessian.has_inf()) {
      Rcpp::Rcout << "The Newton-Raphson algorithm encountered values that "
                  << "produced illegal derivatives." << std::endl;
      happy = false; 
      return lik;
    }
    
    ++iteration;
    vec step = solve(hessian, score.t());
    theta -= step;
    
    double directional_derivative = dot(score, step);
    lik = inference(xx, y, theta, score, hessian);
    lcrit = oldlik - lik;
    
    // if lik decrease step halving 
    if (BAD(lcrit, 0.5e-5)) { 
      if (std::isfinite(lik)) {
        if (directional_derivative < 0) {
          if (fabs(directional_derivative) < 1e-5) return lik;
        }
      }
      
      ++ total_step_halving;
      vec oldtheta = theta + step;
      double step_scale_factor = 1.0;
      while (BAD(lcrit, 0.5e-5) && (step_halving++ <= max_step_halving)) {
        step_scale_factor /= 2.0;
        step *= step_scale_factor;  // halve step size
        theta = oldtheta - step;
        lik = inference(xx, y, theta, score, hessian);
        lcrit = oldlik - lik;
      }
      
      if (rcond(hessian) == 0){ // hessian is singular 
        Rcpp::Rcout << "The Hessian matrix is not positive definite in "
                    << "newton_raphson_min." << std::endl;
        happy = false;
        return lik;
      }
    }
    oldlik = lik;
    if ((step_halving > max_step_halving) || 
        (total_step_halving > max_total_step_halving)) {
      happy = false;
      return lik;
    }
  }
  
  if(iteration == max_iterations || !happy) Rcpp::warning("Maximum likelihood esimation not converged");
  return lik;
}

// sub model 
double Newton(const mat& xx, const vec& y, vec& theta, const uvec& v_k,
              double inference(const mat&, const vec&, const vec&,  const arma::uvec&, rowvec&, mat&),
              int& max_iterations, bool& happy, int marginal = 1){
      
      //Initialize paramters 
      int p = sum(v_k); 
      rowvec score(p);
      mat hessian(p, p);
      theta = theta % v_k;
  
     // Prepre for mle
     int iteration = 0;
     double lcrit = 1, lik = 0;
     double oldlik = inference(xx, y, theta, v_k, score, hessian);
  
     int step_halving = 0, total_step_halving = 0;
     const int max_step_halving = 10, max_total_step_halving = 50;
  
     // Newton update
    while(keep_going(lcrit, 1e-5, iteration, max_iterations)){
       
      if ( score.has_inf() || hessian.has_inf()) {
         Rcpp::Rcout << "The Newton-Raphson algorithm encountered values that "
                     << "produced illegal derivatives." << std::endl;
         happy = false; 
         return lik;
      }
    
      ++iteration;
      vec step = solve(hessian, score.t());
      theta.elem(find(v_k)) -= step;
    
      double directional_derivative = dot(score, step);
      lik = inference(xx, y, theta, v_k, score, hessian);
      lcrit = oldlik - lik;
    
    // if lik decrease step halving 
      if (BAD(lcrit, 0.5e-5)) { 
        if (std::isfinite(lik)) {
          if (directional_derivative < 0) {
            if (fabs(directional_derivative) < 1e-5) return lik;
          }
        }
      
      ++ total_step_halving;
      vec oldtheta = theta.elem(find(v_k)) + step;
      double step_scale_factor = 1.0;
      while (BAD(lcrit, 0.5e-5) && (step_halving++ <= max_step_halving)) {
          step_scale_factor /= 2.0;
          step *= step_scale_factor;  // halve step size
          theta.elem(find(v_k)) = oldtheta - step;
          lik = inference(xx, y, theta, v_k, score, hessian);
          lcrit = oldlik - lik;
        }
      
      if (rcond(hessian) == 0){ // hessian is singular 
        Rcpp::Rcout << "The Hessian matrix is not positive definite in "
                    << "newton_raphson_min." << std::endl;
        happy = false;
        return lik;
      }
    }
      oldlik = lik;
    if ((step_halving > max_step_halving) || 
        (total_step_halving > max_total_step_halving)) {
      happy = false;
      return lik;
    }
  }
    
    
  
  if(iteration == max_iterations || !happy) Rcpp::warning("Maximum likelihood esimation not converged"); 
  return lik;
  }

// full model
Rcpp::List Regression(const arma::mat& xx, const arma::vec& y, int marginal, int maxit){
  int p;
  if(marginal == 1){ // Poisson
    p = xx.n_cols;
  }
  if(marginal == 2){// Negative Binomial
    p = xx.n_cols + 1;
  }
  vec theta(p); theta.zeros();
  theta(0) = log(sum(y)/y.size());
  bool happy = true; 
  double l;
  if(marginal == 1)
    l = Newton(xx, y, theta, pois_infc, maxit, happy); // override maxit to iterations left

  if(marginal == 2){
    theta(p - 1) = 1;
    l = Newton(xx, y, theta, nbinom_infc, maxit, happy, 2);
  }
 
  return Rcpp::List::create(Rcpp::Named("likelihood") = l, 
                            Rcpp::Named("paramters") = theta,
                            Rcpp::Named("happyending") = happy);
}

//sub model 
Rcpp::List Regression(const arma::mat& xx, const arma::vec& y, arma::vec theta, const arma::uvec& v_k, int marginal, int maxit){
  
  bool happy = true; 
  double l;
  if(marginal == 1){
    l = Newton(xx, y, theta, v_k, pois_infc, maxit, happy); 
  }else{
    l = Newton(xx, y, theta, v_k, nbinom_infc, maxit, happy, 2);
  }
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = l, 
                            Rcpp::Named("paramters") = theta,
                            Rcpp::Named("happyending") = happy);
}