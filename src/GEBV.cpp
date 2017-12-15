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

// Retrieves GEBVs for RRDsol
// [[Rcpp::export]]
arma::mat gebvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop){
  arma::mat a = RRsol.slot("markerEff");
  arma::mat d = RRsol.slot("domEff");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 RRsol.slot("lociPerChr"),
                 RRsol.slot("lociLoc"));
  arma::mat output = geno*a+getDomGeno(geno)*d;
  return output;
}

// [[Rcpp::export]]
arma::mat gebvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  bool female, bool isSCAsol=false){
  arma::mat a;
  if(female){
    a = Rcpp::as<arma::mat>(sol.slot("femaleEff"));
  }else{
    a = Rcpp::as<arma::mat>(sol.slot("maleEff"));
  }
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 sol.slot("lociPerChr"),
                 sol.slot("lociLoc"));
  arma::mat X = arma::conv_to<arma::mat>::from(geno);
  if(isSCAsol) X -= 1;
  return X*a;
}

// [[Rcpp::export]]
arma::mat gebvSCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  bool isSCAsol=true){
  arma::mat a1 = sol.slot("femaleEff");
  arma::mat a2 = sol.slot("maleEff");
  arma::mat a3;
  if(isSCAsol){
    a3 = Rcpp::as<arma::mat>(sol.slot("scaEff"));
  }
  arma::Mat<unsigned char> geno;
  geno = getOneHaplo(pop.slot("geno"), 
                     sol.slot("lociPerChr"),
                     sol.slot("lociLoc"), 
                     1);
  arma::mat X1 = arma::conv_to<arma::mat>::from(geno);
  X1 = X1*2;
  if(isSCAsol) X1 -= 1;
  geno = getOneHaplo(pop.slot("geno"), 
                     sol.slot("lociPerChr"),
                     sol.slot("lociLoc"), 
                     2);
  arma::mat X2 = arma::conv_to<arma::mat>::from(geno);
  X2 = X2*2;
  if(isSCAsol) X2 -= 1;
  arma::mat output;
  if(isSCAsol){
    output = X1*a1 + X2*a2 + (X1%X2)*a3;
  }else{
    output = X1*a1 + X2*a2;
  }
  return output;
}
