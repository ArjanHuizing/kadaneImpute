#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can

// [[Rcpp::export]]
vec RIEPS(mat x, mat sigma, int its = 1, double crit = 0.01) {
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
  // start interating
  vec resids = zeros<vec>(m);
  vec old_resids = zeros<vec>(m);
  for(j = 0; j < its; j++){
    // predict obs
    resids = zeros<vec>(m); // reset resids
    for(i = 0; i < m; i++){ // recalc means
      col = x_imp.col(i);
      means[i] = mean(col);
    }
    sigma = cov(x_imp); // recalc sigma
               
    for(i = 0; i < n; i++){
      row = conv_to< vec >::from(x.row(i));
      row_imp = conv_to< vec >::from(x_imp.row(i));
      miss = find_nonfinite(row);
      avail = find_finite(row);
      if(m > miss.n_elem){
        row_imp(avail) = means(avail) + (sigma(avail, miss) * 
          inv_sympd(sigma(miss, miss))) * (row_imp(miss) - means(miss));
        resids(avail) += pow(row(avail) - row_imp(avail).t(), 2) / (nmiss(avail) - 1);  
      }
    }
    // compare
    if((max(resids - old_resids)) > crit){
      if(j < (its-1)){old_resids = resids;}
      
      for(i = 0; i < n; i++){
        row = conv_to< vec >::from(x.row(i));
        row_imp = conv_to< vec >::from(x_imp.row(i));
        miss = find_nonfinite(row);
        if(m > miss.n_elem){
          row_imp(miss) += (randn<vec>(miss.n_elem) * resids(miss));
          x_imp.row(i) = row_imp.t();
        }
      }
    } else {
      break;
    }
  }

  // return
  return resids;
}


/*** R
set.seed(3519)
RIEPS(data, sigma, its = 10)

*/

