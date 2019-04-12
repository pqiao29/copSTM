#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

vec rep_each(const vec& x, const int each) {
  std::size_t n = x.n_elem;
  std::size_t n_out = n*each;
  vec res(n_out);
  auto begin = res.begin();
  for (std::size_t i = 0, ind = 0; i < n; ind += each, ++i) {
    auto start = begin + ind;
    auto end = start + each;
    std::fill(start, end, x[i]);
  }
  return res;
}
mat mat_vec_same_len(mat mt1, vec v1){
  int t = 0;
  for(int i = 0; i != mt1.n_cols; ++i){
    for(int j = 0; j != mt1.n_rows; ++j){
      mt1(j,i) = mt1(j,i) * v1(t);
      ++t;
    }
  }
  return(mt1);
}
vec pmax_c(double a, vec b){
  vec c(b.n_elem);
  for(int i = 0; i != b.n_elem; ++i){
    c(i) = std::max(a, b(i));
  }
  return c;
}

// [[Rcpp::plugins(cpp11)]]
void nearPD_cpp(arma::mat& X // require X: square && symmetry
                  , bool corr = false
                  , double eig_tol   = 1e-6 // defines relative positiveness of eigenvalues compared to largest
                  , double conv_tol  = 1e-7 // convergence tolerance for algorithm
                  , double posd_tol  = 1e-8 // tolerance for enforcing positive definiteness
                  , int maxit = 200 // maximum number of iterations allowed
){
  int n = X.n_cols;
  vec diagX0 = X.diag();
  
  mat D_S; D_S.zeros(n, n);
  
  int iter = 0 ;
  bool converged = false; 
  double conv = R_PosInf;
  mat Y, R, Q;
  vec d;
  
  while (iter < maxit && !converged) {
    Y = X;
    R = Y - D_S;
    
    eig_sym(d, Q, R);

    uvec p = (d > eig_tol * d(n - 1));
    if(sum(p)==0){
      Rcout << "Matrix seems negative semi-definite" << std::endl;
      return;
    } 
    
    uvec p_indexes(sum(p));
    int p_i_i = 0;
    for(int i = 0;i != p.n_elem; ++i){
      if(p(i)){
        p_indexes(p_i_i)=i;
        ++p_i_i;
      }
    }
    
    Q = Q.cols(p_indexes);
    X = mat_vec_same_len(Q,rep_each(d.elem(p_indexes),Q.n_rows))*Q.t();
    D_S = X - R;
    
    if(corr){
      X.diag().ones(); //set diagnols as ones
    } 
    else {
      X.diag() = diagX0;
    } 
    conv = norm(Y-X,"inf")/norm(Y,"inf");
    ++iter;
    converged = (conv <= conv_tol);
  }
  
  if(!converged){ 
    Rcpp::Rcout << "did not converge! " <<std::endl;
  }
  
  eig_sym(d, Q, X);

  double Eps = posd_tol * std::abs(d(n - 1));
  
  if (d(0) < Eps){
    for(int i = 0; i != n; ++i){
      if(d(i) < Eps){
        d(i) = Eps;
      }
    }

    vec o_diag = X.diag();
    
    mat Q_t = Q.t();
    for(int i = 0; i != Q_t.n_cols; ++i){
      Q_t.col(i) = Q_t.col(i) % d;
    }
    X = Q * Q_t;
    vec D = sqrt(pmax_c(Eps, o_diag)/X.diag());
    
    for(int j = 0; j != X.n_cols; ++j){
      X.col(j) = X.col(j) % D;
    }
    for(int i = 0; i != X.n_rows; ++i){
      X.row(i) = X.row(i) % conv_to<rowvec>::from(D);
    }
    
  }
  
  if(corr) {
    X.diag().ones(); //set diag as ones
  }
  else {
    X.diag() = diagX0;
  } 
  
}