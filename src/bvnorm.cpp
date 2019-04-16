#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <vector>

double Phi2diag( const double& x,
                 const double& a,  // 1 - rho
                 const double& px, // Phi( x )
                 const double& pxs) // Phi( lambda( rho ) * x )
{
  if( a <= 0.0 ) return px;  // rho == 1
  if( a >= 1.0 ) return px * px;  // rho == 0
  
  double b = 2.0 - a, sqrt_ab = sqrt( a * b );
  double asr = ( a > 0.1 ? asin( 1.0 - a ) : acos( sqrt_ab ) );
  double comp = px * pxs;
  if( comp * ( 1.0 - a - 6.36619772367581343e-001 * asr ) < 5e-17 )
    return b * comp;
  
  double tmp = 1.25331413731550025 * x;
  double a_coeff = a * x * x / b;
  double a_even = -tmp * a;
  double a_odd = -sqrt_ab * a_coeff;
  double b_coeff = x * x;
  double b_even = tmp * sqrt_ab;
  double b_odd = sqrt_ab * b_coeff;
  double d_coeff = 2.0 * x * x / b;
  double d_even = ( 1.0 - a ) * 1.57079632679489662 - asr;
  double d_odd = tmp * ( sqrt_ab - a );
  
  double res = 0.0, res_new = d_even + d_odd;
  int k = 2;
  while( res != res_new )
  {
    d_even = ( a_odd + b_odd + d_coeff * d_even ) / k;
    a_even *= a_coeff / k;
    b_even *= b_coeff / k;
    ++k;
    a_odd *= a_coeff / k;
    b_odd *= b_coeff / k;
    d_odd = ( a_even + b_even + d_coeff * d_odd ) / k;
    ++k;
    res = res_new;
    res_new += d_even + d_odd;
  }
  res *= exp( -x * x / b ) * 1.591549430918953358e-001;
  return std::max( ( 1.0 + 6.36619772367581343e-001 * asr ) * comp,
              b * comp - std::max( 0.0, res ) );
}

double Phi(double value)  // pnorm
{
  return 0.5 * erfc(-value * M_SQRT1_2);
}

double Phi2help( const double& x,
                 const double& y,
                 const double& rho ){
  
  if( x == 0.0 ) return ( y >= 0.0 ? 0.0 : 0.5 );
  
  double s = sqrt( ( 1.0 - rho ) * ( 1.0 + rho ) );
  double a = 0.0, b1 = -fabs( x ), b2 = 0.0;
  if( rho > 0.99 )
  {
    double tmp = sqrt( ( 1.0 - rho ) / ( 1.0 + rho ) );
    b2 = -fabs( ( x - y ) / s - x * tmp );
    a = (( x - y ) / x / s - tmp) * (( x - y ) / x / s - tmp);
  } else if( rho < -0.99 )
  {
    double tmp = sqrt( ( 1.0 + rho ) / ( 1.0 - rho ) );
    b2 = -fabs( ( x + y ) / s - x * tmp );
    a = ( ( x + y ) / x / s - tmp ) * ( ( x + y ) / x / s - tmp );
  }else {
    b2 = -fabs( rho * x - y ) / s;
    a = ( b2 / x ) * ( b2 / x );
  }
  
  double p1 = Phi( b1 ), p2 = Phi( b2 ); // cum. standard normal
  double q = 0.0;
  if( a <= 1.0 )
  {q = 0.5 * Phi2diag( b1, 2.0 * a / ( 1.0 + a ), p1, p2 );}
  else
  {q = p1 * p2 - 0.5 * Phi2diag( b2, 2.0 / ( 1.0 + a ), p2, p1 );}
  
  int c1 = ( y / x >= rho );
  int c2 = ( x < 0.0 );
  int c3 = c2 && ( y >= 0.0 );
  
  return ( c1 && c3 ? q - 0.5 : c1 && c2 ? q
             : c1 ? 0.5 - p1 + q
             : c3 ? p1 - q - 0.5
             : c2 ? p1 - q
             : 0.5 - q );
  
}


double dbnorm_inf( const double& x,
             const double& y,
             const double& rho ){
    
    if(std::isinf(x) || std::isinf(y))
        return 0;
    
  if( ( 1.0 - rho ) * ( 1.0 + rho ) <= 0.0 ){
    if( rho > 0.0 ){
      return Phi( std::min( x, y ) );
    } else {
      return std::max( 0.0, std::min( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ); }
  }

  if( x == 0.0 && y == 0.0 ){ if( rho > 0.0 ){
    return Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 );
  } else {
    return 0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 );
  }
  }

  return max( 0.0,
              min( 1.0,
                   Phi2help( x, y, rho ) + Phi2help( y, x, rho ) ) );
}



double dbnorm( const std::vector<double>& x,
               const std::vector<double>& y,
               const double& rho ){
  
  /* requirement: x.size() == 2 && y.size() == 2
  x[0] <= x[1] && y[0] <= y[1]
  */
  
  if(x[0] == x[1] || y[0] == y[1]){
    return 0;
  }else{
    if(std::isinf(x[1]) && std::isinf(y[1])){
      return 1;
    }else if(std::isinf(x[1])){
      return Phi(y[1]) - Phi(y[0]);
    }else if(std::isinf(y[1])){
      return Phi(x[1]) - Phi(x[0]);
    }else{
      return dbnorm_inf(x[1], y[1], rho) + dbnorm_inf(x[0], y[0], rho) - dbnorm_inf(x[0], y[1], rho) - dbnorm_inf(x[1], y[0], rho);
    }
  }
  
}