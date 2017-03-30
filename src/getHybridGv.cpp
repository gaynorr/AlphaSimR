// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Retrieves hybrid geno under the assumption lines are completely inbred
arma::Mat<unsigned char> getHybridGeno(const Rcpp::S4& trait, 
                                       const Rcpp::S4& fPop, arma::uvec fPar,
                                       const Rcpp::S4& mPop, arma::uvec mPar){
  // R to C++
  fPar = fPar-1;
  mPar = mPar-1;
  arma::Mat<unsigned char> geno;
  geno = (getGeno(fPop.slot("geno"), 
                  trait.slot("lociPerLoc"),
                  trait.slot("lociLoc"))).rows(fPar);
  geno += (getGeno(mPop.slot("geno"), 
                   trait.slot("lociPerLoc"),
                   trait.slot("lociLoc"))).rows(mPar);
  geno = geno/2;
  return geno;
}

// Retrieves genetic values for TraitA
// [[Rcpp::export]]
arma::vec getHybridGvA(const Rcpp::S4& trait, 
                       const Rcpp::S4& fPop, arma::uvec& fPar,
                       const Rcpp::S4& mPop, arma::uvec& mPar){
  
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno = 
    getHybridGeno(trait,fPop,fPar,mPop,mPar);
  return calcGvA(geno, a, intercept);
}

// Retrieves genetic values for TraitAD
// [[Rcpp::export]]
arma::vec getHybridGvAD(const Rcpp::S4& trait, 
                        const Rcpp::S4& fPop, arma::uvec& fPar,
                        const Rcpp::S4& mPop, arma::uvec& mPar){
  arma::vec a = trait.slot("addEff");
  arma::vec d = trait.slot("domEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno = 
    getHybridGeno(trait,fPop,fPar,mPop,mPar);
  return calcGvAD(geno, a, d, intercept);
}
