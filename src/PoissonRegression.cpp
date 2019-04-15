#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 


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

// require: positive y, xx.n_rows == y.size(), xx.n_cols == theta.size()
double pois_inference(const mat& xx, const vec& y, const vec& theta
                        , rowvec& score, mat& hessian) 
{
  vec lmd = exp(xx * theta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    lik -= R::dpois(y(i), lmd(i), true);
    rowvec tmp_row = xx.row(i);
    score -= (y(i) - lmd(i)) * tmp_row;
    hessian += lmd(i) * (tmp_row.t() * tmp_row);
  }
  return lik;
}

double pois_inference(const mat& xx, const vec& y, const vec& theta, 
                      const uvec& v, rowvec& score, mat& hessian) 
{
  vec lmd = exp(xx * theta);
  double lik = 0;
  score.zeros(); hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    lik -= R::dpois(y(i), lmd(i), true);
    vec tmp_row = xx.row(i).t();  
    tmp_row = tmp_row.elem(find(v));
    score -= (y(i) - lmd(i)) * tmp_row.t();
    hessian += lmd(i) * (tmp_row * tmp_row.t());
  }
  return lik;
}

void get_hessian(const mat& xx, const vec& y, const vec& theta, mat& hessian){
  vec lmd = exp(xx * theta);
  hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    rowvec tmp_row = xx.row(i);
    hessian += lmd(i) * (tmp_row.t() * tmp_row);
  }
}

void get_hessian(const mat& xx, const vec& y, const vec& theta, const uvec& v, mat& hessian){
  vec lmd = exp(xx * theta);
  hessian.zeros();
  for(int i = 0; i != y.size(); ++i){
    vec tmp_row = xx.row(i).t();  
    tmp_row = tmp_row.elem(find(v));
    hessian += lmd(i) * (tmp_row * tmp_row.t());
  }
}

double Poisson_Newton(const mat& xx, const vec& y, vec& theta, 
                      int& max_iterations, bool& happy){
  //Initialize paramters 
  int p = theta.size();
  rowvec score(p);
  mat hessian(p, p);
  
  // Prepre for mle
  int iteration = 0;
  double lcrit = 1, lik = 0;
  double oldlik = pois_inference(xx, y, theta, score, hessian);
  
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
          lik = pois_inference(xx, y, theta, score, hessian);
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
        lik = pois_inference(xx, y, theta, score, hessian);
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
  
  if(iteration == max_iterations || !happy) Rcpp::warning("Maximum likelihood esimation not converged"); else max_iterations -= iteration;
  
  return lik;
  }

double Poisson_Newton(const mat& xx, const vec& y, const vec& theta, const uvec& v_k, 
                      const int max_iterations, bool happy = true){
  //Initialize paramters 
  int p = sum(v_k);
  rowvec score(p);
  mat hessian(p, p);
  vec tmp_theta = theta % v_k;
  
  // Prepre for mle
  int iteration = 0;
  double lcrit = 1, lik = 0;
  double oldlik = pois_inference(xx, y, tmp_theta, v_k, score, hessian);
  
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
    tmp_theta.elem(find(v_k)) -= step;
    
    double directional_derivative = dot(score, step);
    lik = pois_inference(xx, y, tmp_theta, v_k, score, hessian);
    lcrit = oldlik - lik;
    
    
    // if lik decrease step halving 
    if (BAD(lcrit, 0.5e-5)) { 
      if (std::isfinite(lik)) {
        if (directional_derivative < 0) {
          if (fabs(directional_derivative) < 1e-5) return lik;
        }
      }
      
      ++ total_step_halving;
      vec oldtheta = tmp_theta.elem(find(v_k)) + step;
      double step_scale_factor = 1.0;
      while (BAD(lcrit, 0.5e-5) && (step_halving++ <= max_step_halving)) {
        step_scale_factor /= 2.0;
        step *= step_scale_factor;  // halve step size
        tmp_theta.elem(find(v_k)) = oldtheta - step;
        lik = pois_inference(xx, y, tmp_theta, v_k, score, hessian);
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

double Poisson_Newton(const mat& xx, const vec& y, vec& theta, const uvec& v_k, 
                    const int max_iterations, bool happy = true){
  //Initialize paramters 
  int p = sum(v_k);
  rowvec score(p);
  mat hessian(p, p);
  theta = theta % v_k;
  
  // Prepre for mle
  int iteration = 0;
  double lcrit = 1, lik = 0;
  double oldlik = pois_inference(xx, y, theta, v_k, score, hessian);
  
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
    lik = pois_inference(xx, y, theta, v_k, score, hessian);
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
        lik = pois_inference(xx, y, theta, v_k, score, hessian);
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
    }
  }
  
  // No warnings generated, since this function is only called for selected model
  
  return lik;
}


Rcpp::List PoissonRegression(const arma::mat& xx, const arma::vec& y, int& maxit){
  int p = xx.n_cols;
  vec theta(p); theta.zeros();
      theta(0) = sum(y)/y.size();
  bool happy = true; 
  double l = Poisson_Newton(xx, y, theta, maxit, happy); // override maxit to iterations left
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = l, 
                            Rcpp::Named("paramters") = theta, 
                            Rcpp::Named("happyending") = happy
                             );
}
