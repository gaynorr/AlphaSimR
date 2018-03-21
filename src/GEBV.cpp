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
  return arma::conv_to<arma::mat>::from(geno)*a;
}

// Retrieves GEGVs for RRDsol
// [[Rcpp::export]]
arma::mat gegvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop){
  arma::mat a = RRsol.slot("addEff");
  arma::mat d = RRsol.slot("domEff");
  double b = RRsol.slot("hetCov");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 RRsol.slot("lociPerChr"),
                 RRsol.slot("lociLoc"));
  arma::Mat<unsigned char> genoD = getDomGeno(geno);
  arma::mat het = mean(arma::conv_to<arma::mat>::from(genoD),1);
  return arma::conv_to<arma::mat>::from(geno)*a+
    arma::conv_to<arma::mat>::from(genoD)*d+het*b;
}

// Retrieves GEBVs for RRDsol using population specific p
// [[Rcpp::export]]
arma::mat gebvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop){
  arma::mat a = RRsol.slot("addEff");
  arma::mat d = RRsol.slot("domEff");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 RRsol.slot("lociPerChr"),
                 RRsol.slot("lociLoc"));
  arma::mat p = arma::mean(arma::conv_to<arma::mat>::from(geno),0)/2;
  return arma::conv_to<arma::mat>::from(geno)*(a+d%(1-2*p.t()));
}

// [[Rcpp::export]]
arma::mat gebvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  bool female){
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
  return arma::conv_to<arma::mat>::from(geno)*a;
}

// [[Rcpp::export]]
arma::mat gebvSCA_GCA(const Rcpp::S4& sol, const Rcpp::S4& pop){
  arma::mat a1 = sol.slot("femaleEff");
  arma::mat a2 = sol.slot("maleEff");
  arma::Mat<unsigned char> geno1,geno2;
  geno1 = getOneHaplo(pop.slot("geno"), 
                      sol.slot("lociPerChr"),
                      sol.slot("lociLoc"), 
                      1);
  geno2 = getOneHaplo(pop.slot("geno"), 
                      sol.slot("lociPerChr"),
                      sol.slot("lociLoc"), 
                      2);
  return 2*(arma::conv_to<arma::mat>::from(geno1)*a1+
            arma::conv_to<arma::mat>::from(geno2)*a2);
}

// [[Rcpp::export]]
arma::mat gebvSCA_SCA(const Rcpp::S4& sol, const Rcpp::S4& pop){
  arma::mat a1 = sol.slot("a1");
  arma::mat a2 = sol.slot("a2");
  arma::mat d = sol.slot("d");
  double b = sol.slot("hetCov");
  arma::Mat<unsigned char> geno1,geno2,genoD;
  geno1 = getOneHaplo(pop.slot("geno"), 
                      sol.slot("lociPerChr"),
                      sol.slot("lociLoc"), 
                      1);
  geno2 = getOneHaplo(pop.slot("geno"), 
                      sol.slot("lociPerChr"),
                      sol.slot("lociLoc"), 
                      2);
  genoD = getDomGeno(geno1+geno2);
  arma::mat het = mean(arma::conv_to<arma::mat>::from(genoD),1);
  return 2*(arma::conv_to<arma::mat>::from(geno1)*a1+
            arma::conv_to<arma::mat>::from(geno1)*a1)+
            arma::conv_to<arma::mat>::from(genoD)*d+het*b;
}
