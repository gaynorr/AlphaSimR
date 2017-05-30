// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"


// Retrieves GEBVs for RRsol
// [[Rcpp::export]]
arma::mat gebvRR(const Rcpp::S4& RRsol, const Rcpp::S4& pop){
  arma::mat a = RRsol.slot("markerEff");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 RRsol.slot("lociPerChr"),
                 RRsol.slot("lociLoc"));
  arma::mat output = geno*a;
  return output;
}

// [[Rcpp::export]]
arma::mat gebvGCA(const Rcpp::S4& GCAsol, const Rcpp::S4& pop, 
                  bool female){
  arma::mat a;
  if(female){
    a = Rcpp::as<arma::mat>(GCAsol.slot("femaleEff"));
  }else{
    a = Rcpp::as<arma::mat>(GCAsol.slot("maleEff"));
  }
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 GCAsol.slot("lociPerChr"),
                 GCAsol.slot("lociLoc"));
  arma::mat X = arma::conv_to<arma::mat>::from(geno);
  X = X/2;
  return X*a;
}

// [[Rcpp::export]]
arma::mat gebvSCA(const Rcpp::S4& SCAsol, const Rcpp::S4& pop){
  arma::mat a1 = SCAsol.slot("femaleEff");
  arma::mat a2 = SCAsol.slot("maleEff");
  arma::mat a3 = SCAsol.slot("scaEff");
  arma::Mat<unsigned char> geno;
  geno = getOneHaplo(pop.slot("geno"), 
                     SCAsol.slot("lociPerChr"),
                     SCAsol.slot("lociLoc"), 
                     1);
  arma::mat X1 = arma::conv_to<arma::mat>::from(geno);
  X1 = X1*2-1;
  geno = getOneHaplo(pop.slot("geno"), 
                     SCAsol.slot("lociPerChr"),
                     SCAsol.slot("lociLoc"), 
                     2);
  arma::mat X2 = arma::conv_to<arma::mat>::from(geno);
  X2 = X2*2-1;
  arma::mat output = X1*a1 + X2*a2 + (X1%X2)*a3;
  return output;
}
