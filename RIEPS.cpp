#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

// The RIEPS algorithm (Rassler, 2002) creates a second set of predictions, predicting observed
// values using imputed ones. These are then used to estimate the residual
// variance.

// The initial estimate is biased, and underestimates the amount of variance.
// This behaviour is adressed by iterating the second set of predictions with
// updated residual variances (Kiesl & Rassler, 2009).

// [[Rcpp::export]]
mat RIEPS(mat x, mat sigma, int its = 100, double crit = 0.005) {
  int n = x.n_rows;
  double nrow = x.n_rows;
  int m = x.n_cols;
  int i, j;
  vec col(n);
  vec row(m);
  vec row_imp(m);
  vec means(m);
  vec nmiss(m);
  uvec avail(n);
  uvec miss(n);
  mat x_imp = x;
  
  // means
  for(j = 0; j < m; j++){
    col = x.col(j);
    miss = find_nonfinite(col);
    avail = find_finite(col);
    means[j] = mean(col.elem(avail));
  }
  // predict miss
  for(i = 0; i < n; i++){
    row = conv_to< vec >::from(x.row(i));
    row_imp = conv_to< vec >::from(x.row(i));
    miss = find_nonfinite(row);
    avail = find_finite(row);
    if(m > avail.n_elem){
      row_imp(miss) = means(miss) + (sigma(miss, avail) * 
        inv_sympd(sigma(avail, avail))) * (row(avail) - means(avail));
      x_imp.row(i) = row_imp.t();  
      nmiss(miss) += 1;
    }
  }
  
  // start iterating
  vec resids = zeros<vec>(m);
  vec old_resids = zeros<vec>(m);
  mat x_s = x_imp;
  for(j = 0; j < its; j++){
    // predict obs
    resids = zeros<vec>(m); // reset resids
    for(i = 0; i < m; i++){ // recalc means
      col = x_s.col(i);
      means[i] = mean(col);
    }
    sigma = cov(x_s); // recalc sigma
               
    for(i = 0; i < n; i++){
      row = conv_to< vec >::from(x.row(i));
      row_imp = conv_to< vec >::from(x_s.row(i));
      miss = find_nonfinite(row);
      avail = find_finite(row);
      if(m > miss.n_elem){
        row_imp(avail) = means(avail) + (sigma(avail, miss) * 
          inv_sympd(sigma(miss, miss))) * (row_imp(miss) - means(miss));
        
        resids(avail) += pow((row(avail) - row_imp(avail)), 2) / (nmiss(avail) - 1);  
      }
    }
    // compare
    if((max(resids/old_resids)) > (1 + crit)){
      old_resids = resids;
      // add resids
      for(i = 0; i < n; i++){
        row = conv_to< vec >::from(x.row(i));
        row_imp = conv_to< vec >::from(x_imp.row(i));
        miss = find_nonfinite(row);
        if(m > miss.n_elem){
          row_imp(miss) += (sqrt(resids(miss)) % randn<vec>(miss.n_elem));
          x_s.row(i) = row_imp.t();
        }
      }
    } else {
      break;
    }
  }
  // return
  return x_s;
}