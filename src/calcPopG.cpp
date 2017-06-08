// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// [[Rcpp::export(.calcPopG)]]
arma::fmat calcPopG(const arma::field<arma::Cube<unsigned char> >& geno, 
                                  const arma::ivec& lociPerChr,
                                  arma::uvec lociLoc){
  arma::fmat X = arma::conv_to<arma::fmat>::from(getGeno(geno, lociPerChr, lociLoc));
  arma::frowvec p = mean(X,0)/2.0;
  X.each_row() -= 2*p;
  arma::fmat G = X*X.t();
  G = G/(2.0*sum(p%(1-p)));
  return G;
}

// [[Rcpp::export(.calcPopGIbs)]]
arma::fmat calcPopGIbs(const arma::field<arma::Cube<unsigned char> >& geno, 
                       const arma::ivec& lociPerChr,
                       arma::uvec lociLoc){
  arma::fmat X = arma::conv_to<arma::fmat>::from(getGeno(geno, lociPerChr, lociLoc));
  X -= 1.0;
  arma::fmat G = (X*X.t())/X.n_cols + 1.0;
  return G;
  return G;
}
