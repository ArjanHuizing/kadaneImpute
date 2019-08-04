#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

// Estimate the covariance matrix for a dataset with missing values using
// maximum likelihood estimation. Values are obtained through a
// Expectation-Maximization algorithm.

// [[Rcpp::export]]
mat covMiss(mat x, int its = 100) {
  int n = x.n_rows;
  double nrow = x.n_rows;
  int m = x.n_cols;
  int i, j, k;
  vec means(m);
  vec col(n);
  vec row(m);
  vec row_imp(m);
  uvec avail(n);
  uvec miss(n);
  mat x_imp = x;
  
  // impute missings with means as starting values.
  for(j = 0; j < m; j++){
    col = x.col(j);
    miss = find_nonfinite(col);
    avail = find_finite(col);
    // obtain means and impute
    means[j] = mean(col.elem(avail));
    col.elem(miss).fill(means[j]);
    // assign imputed column
    x_imp.col(j) = col;
    means[j] = mean(col);
  }
  // starting values
  mat sigma = cov(x_imp) * (nrow-1)/nrow;
  mat bias(m, m);  
  
  // EM algorithm
  for(j = 0; j < its; j++){
    bias = zeros<mat>(m, m);
    for(i = 0; i < n; i++){
      row = conv_to< vec >::from(x.row(i));
      row_imp = conv_to< vec >::from(x.row(i));
      miss = find_nonfinite(row);
      avail = find_finite(row);
      if(m > avail.n_elem){
        bias(miss, miss) = bias(miss, miss) + sigma(miss, miss); - 
          sigma(miss, avail) * inv_sympd(sigma(avail, avail)) * sigma(avail, miss);
        
        row_imp(miss) = means(miss) + (sigma(miss, avail) * 
          inv_sympd(sigma(avail, avail))) * (row(avail) - means(avail));
        x_imp.row(i) = row_imp.t();  
      }
    }
    
    // update
    for(k = 0; k < m; k++){
      col = x_imp.col(k);
      means[k] = mean(col);
    }
    sigma = (cov(x_imp) * (nrow-1)/nrow) + (bias/nrow);
  }
  // return
  return sigma;
}