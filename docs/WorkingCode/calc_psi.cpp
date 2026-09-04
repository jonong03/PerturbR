// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat calc_psi_cpp(const arma::mat& R, double n = 1.0) {
  
  int p = R.n_rows;
  
  if (R.n_cols != p)
    stop("R must be a square matrix.");
  
  int k = p * (p - 1) / 2;
  
  // store correlation-index pairs
  arma::imat idx(k, 2);
  
  int z = 0;
  for (int e = 0; e < p - 1; ++e) {
    for (int f = e + 1; f < p; ++f) {
      idx(z, 0) = e;
      idx(z, 1) = f;
      ++z;
    }
  }
  
  arma::mat psi(k, k, arma::fill::zeros);
  
  for (int i = 0; i < k; ++i) {
    
    int e = idx(i, 0);
    int f = idx(i, 1);
    
    double ref = R(e, f);
    
    psi(i, i) =
      (1.0 - ref * ref) *
      (1.0 - ref * ref);
    
    for (int j = 0; j < i; ++j) {
      
      int g = idx(j, 0);
      int h = idx(j, 1);
      
      double tmp1 =
        (R(e,g) - R(e,f) * R(f,g)) *
        (R(f,h) - R(f,g) * R(g,h));
      
      double tmp2 =
        (R(e,h) - R(e,g) * R(g,h)) *
        (R(f,g) - R(f,e) * R(e,g));
      
      double tmp3 =
        (R(e,g) - R(e,h) * R(h,g)) *
        (R(f,h) - R(f,e) * R(e,h));
      
      double tmp4 =
        (R(e,h) - R(e,f) * R(f,h)) *
        (R(f,g) - R(f,h) * R(h,g));
      
      double value =
        0.5 * (tmp1 + tmp2 + tmp3 + tmp4);
      
      psi(i,j) = value;
      psi(j,i) = value;
    }
  }
  
  return psi / n;
}