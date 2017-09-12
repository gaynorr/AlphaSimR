// Functions for adjusting traits to desired initial variance
// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

//Objective function for tuning TraitA
Rcpp::List traitAObj(double tuneValue, Rcpp::List args){
  arma::Mat<unsigned char> geno = args["geno"];
  arma::vec addEff = args["addEff"];
  double varG = args["varG"];
  arma::vec gv = geno*(addEff*tuneValue);
  double intercept = arma::mean(gv);
  double obsVar = arma::var(gv,1); //Population variance
  Rcpp::List output;
  output = Rcpp::List::create(Rcpp::Named("intercept")=intercept);
  return Rcpp::List::create(Rcpp::Named("objective")=fabs(varG-obsVar),
                            Rcpp::Named("output")=output);
}

//Objective function for tuning TraitD
Rcpp::List traitADObj(double tuneValue, Rcpp::List args){
  arma::Mat<unsigned char> geno = args["geno"];
  arma::imat domGeno = args["domGeno"];
  arma::vec addEff = args["addEff"];
  arma::vec domEff = args["domEff"];
  double varG = args["varG"];
  bool useVarA = args["useVarA"];
  arma::vec gv = geno*(addEff*tuneValue)+domGeno*(domEff*tuneValue);
  double intercept = arma::mean(gv);
  double obsVar;
  if(useVarA){
    arma::vec p = (arma::mean(arma::conv_to<arma::mat>::from(geno),0)/2.0).t();
    arma::vec alpha = (addEff*tuneValue)+(domEff*tuneValue)%(1-2*p);
    arma::vec bv = geno*alpha;
    obsVar = arma::var(bv,1);
  }else{
    obsVar = arma::var(gv,1); //Population variance
  }
  Rcpp::List output;
  output = Rcpp::List::create(Rcpp::Named("intercept")=intercept);
  return Rcpp::List::create(Rcpp::Named("objective")=fabs(varG-obsVar),
                            Rcpp::Named("output")=output);
}

//Tunes TraitA for desired varG
// [[Rcpp::export]]
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      double varG){
  return optimize(*traitAObj,
                  Rcpp::List::create(Rcpp::Named("geno")=geno,
                                     Rcpp::Named("addEff")=addEff,
                                     Rcpp::Named("varG")=varG),
                                     1e-10,
                                     1e3);
}

//Tunes TraitAD for desired varG (or varA), tuneTraitA used for inbreds
// [[Rcpp::export]]
Rcpp::List tuneTraitAD(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      arma::vec& domEff,
                      double varG, 
                      bool useVarA){
  return optimize(*traitADObj,
                  Rcpp::List::create(Rcpp::Named("geno")=geno,
                                     Rcpp::Named("domGeno")=getDomGeno(geno),
                                     Rcpp::Named("addEff")=addEff,
                                     Rcpp::Named("domEff")=domEff,
                                     Rcpp::Named("varG")=varG,
                                     Rcpp::Named("useVarA")=useVarA),
                                     1e-10,
                                     1e3);
}





