// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates population variance
// [[Rcpp::export]]
arma::mat popVar(const arma::mat& X) {
  return arma::var(X,1);
}

// Merges geno objects, i.e. fields containing cubes of unsigned char
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > mergeGeno(
    const arma::field<arma::Cube<unsigned char> >& x, 
    const arma::field<arma::Cube<unsigned char> >& y){
  int nChr = x.n_elem;
  arma::field<arma::Cube<unsigned char> > z(nChr);
  for(int i=0; i<nChr; ++i){
    z(i) = arma::join_slices(x(i),y(i));
  }
  return z;
}

/*
 * The functions below are for modifying R objects in place.
 * 
 * The data type of oldValue must match exactly the data type in R or else
 * modifying in place will not occur.
 * 
 * For example, if assignMat is called on an integer matrix the function will
 * not modify in place. This is due to the function casting oldValue from 
 * integer matrix to a numeric matrix prior to assigning newValue.
 */
// [[Rcpp::export]]
void assignMat(arma::mat& oldValue, arma::mat& newValue){
  oldValue = newValue;
  return;
}

// [[Rcpp::export]]
void assignInt(int& oldValue, int& newValue){
  oldValue = newValue;
  return;
}
