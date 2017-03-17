// Functions for adjusting traits to desired initial variance
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "optimize.h"

//Objective function for tuning TraitA
Rcpp::List traitAObj(double tuneValue, Rcpp::List args){
  arma::Mat<unsigned char> geno = args["geno"];
  arma::vec addEff = args["addEff"];
  double varG = args["varG"];
  arma::vec gv = geno*(addEff*tuneValue);
  double intercept = arma::mean(gv);
  double obsVar = arma::var(gv);
  Rcpp::List output;
  output = Rcpp::List::create(Rcpp::Named("intercept")=intercept);
  return Rcpp::List::create(Rcpp::Named("objective")=fabs(varG-obsVar),
                            Rcpp::Named("output")=output);
}

//Tunes TraitA for desired varG
// [[Rcpp::export]]
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      double varG) {
  return optimize(*traitAObj,
                  Rcpp::List::create(Rcpp::Named("geno")=geno,
                                     Rcpp::Named("addEff")=addEff,
                                     Rcpp::Named("varG")=varG),
                                     1e-10,
                                     varG*10);
}

