// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Retrieves hybrid geno under the assumption lines are completely inbred
arma::Mat<unsigned char> getHybridGeno(const Rcpp::S4& trait, 
                                       const Rcpp::S4& motherGeno, arma::uvec mother,
                                       const Rcpp::S4& fatherGeno, arma::uvec father){
  // R to C++
  mother = mother-1;
  father = father-1;
  arma::Mat<unsigned char> geno;
  geno = (getGenoT(motherGeno.slot("geno"), 
                   trait.slot("lociPerChr"),
                   trait.slot("lociLoc"))).cols(mother);
  geno += (getGenoT(fatherGeno.slot("geno"), 
                    trait.slot("lociPerChr"),
                    trait.slot("lociLoc"))).cols(father);
  geno = geno/2;
  return geno;
}

// [[Rcpp::export]]
arma::field<arma::vec> getHybridGv(const Rcpp::S4& trait, 
                                   const Rcpp::S4& motherGeno, arma::uvec& mother,
                                   const Rcpp::S4& fatherGeno, arma::uvec& father, 
                                   int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  if(hasGxe){
    output.set_size(2);
  }else{
    output.set_size(1);
  }
  arma::Mat<unsigned char> geno;
  geno = getHybridGeno(trait,motherGeno,mother,fatherGeno,father);
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  if(hasD){
    arma::vec d = trait.slot("domEff");
    output(0) = calcGvAD(geno, a, d, intercept, nThreads);
  }else{
    output(0) = calcGvA(geno, a, intercept, nThreads);
  }
  if(hasGxe){
    arma::vec g = trait.slot("gxeEff");
    double gxeInt = trait.slot("gxeInt");
    output(1) = calcGvA(geno, g, gxeInt, nThreads);
  }
  return output;
}
