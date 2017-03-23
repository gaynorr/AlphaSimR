// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//Calculates population variance
// [[Rcpp::export]]
arma::mat popVar(arma::mat& X) {
  return arma::var(X,1);
}
