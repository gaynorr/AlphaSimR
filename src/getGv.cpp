// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "getGeno.h"

// [[Rcpp::export]]
arma::vec getGvA(Rcpp::S4& trait, Rcpp::S4& pop){
  arma::vec output;
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop, trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  output = geno*a;
  return output+intercept;
}

// [[Rcpp::export]]
arma::vec getGvAD(Rcpp::S4& trait, Rcpp::S4& pop){
  arma::vec output;
  arma::vec a = trait.slot("addEff");
  arma::vec d = trait.slot("domEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop, trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  output = geno*a+getDomGeno(geno)*d;
  return output+intercept;
}