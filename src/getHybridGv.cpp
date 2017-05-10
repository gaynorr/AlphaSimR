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
                  trait.slot("lociPerChr"),
                  trait.slot("lociLoc"))).rows(fPar);
  geno += (getGeno(mPop.slot("geno"), 
                   trait.slot("lociPerChr"),
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

// Retrieves genetic values for TraitAG
// [[Rcpp::export]]
arma::vec getHybridGvAG(const Rcpp::S4& trait, 
                        const Rcpp::S4& fPop, arma::uvec& fPar,
                        const Rcpp::S4& mPop, arma::uvec& mPar,
                        double z){
  arma::vec a = trait.slot("addEff");
  arma::vec x = trait.slot("gxeEff");
  a = a+x*z;
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

// Retrieves genetic values for TraitADG
// [[Rcpp::export]]
arma::vec getHybridGvADG(const Rcpp::S4& trait, 
                         const Rcpp::S4& fPop, arma::uvec& fPar,
                         const Rcpp::S4& mPop, arma::uvec& mPar,
                         double z){
  arma::vec aOld = trait.slot("addEff");
  arma::vec d = trait.slot("domEff");
  arma::vec x = trait.slot("gxeEff");
  arma::vec a = aOld+x*z;
  d = d%(a/aOld);
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno = 
    getHybridGeno(trait,fPop,fPar,mPop,mPar);
  return calcGvAD(geno, a, d, intercept);
}