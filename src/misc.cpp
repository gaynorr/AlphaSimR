// These function are called by R, but are not listed in the package namespace
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

// Calculates minor allele frequency on a single chromsome
// Requires bi-allelic markers, but works for any ploidy
// [[Rcpp::export]]
arma::vec calcChrMinorFreq(const arma::Cube<unsigned char>& geno,
                           int ploidy){
  arma::Mat<unsigned char> tmp = arma::sum(geno,1);
  arma::vec output = arma::mean(arma::conv_to<arma::mat>::from(tmp),
                                1)/ploidy;
  return 0.5-arma::abs(output-0.5);
}

// [[Rcpp::export]]
arma::imat convToImat(const arma::Mat<unsigned char>& X){
  return arma::conv_to<arma::imat>::from(X);
}

