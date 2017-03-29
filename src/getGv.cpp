// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// [[Rcpp::export]]
arma::vec calcGvA(const arma::Mat<unsigned char>& geno,
                  const arma::vec& a, double intercept){
  arma::vec output = geno*a;
  return output+intercept;
}

// Retrieves genetic value for TraitA
// [[Rcpp::export]]
arma::vec getGvA(Rcpp::S4& trait, Rcpp::S4& pop){
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop, trait);
  return calcGvA(geno, a, intercept);
}

// [[Rcpp::export]]
arma::vec calcGvAD(const arma::Mat<unsigned char>& geno,
                   const arma::vec& a, const arma::vec& d,
                   double intercept){
  arma::vec output = geno*a+getDomGeno(geno)*d;
  return output+intercept;
}

// Retrieves genetic value for TraitAD
// [[Rcpp::export]]
arma::vec getGvAD(Rcpp::S4& trait, Rcpp::S4& pop){
  arma::vec a = trait.slot("addEff");
  arma::vec d = trait.slot("domEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop, trait);
  return calcGvAD(geno,a,d,intercept);
}

// [[Rcpp::export]]
Rcpp::List calcGenParam(Rcpp::S4& trait, Rcpp::S4& pop,
                        arma::vec a, arma::vec d){
  int nInd = pop.slot("nInd");
  arma::vec bv(nInd);
  arma::vec dd(nInd);
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop, trait);
  arma::mat X = arma::conv_to<arma::mat>::from(geno);
  arma::rowvec p = arma::mean(X,0)/2.0;
  arma::vec alpha = a+d%(1-2*p.t()); //allele subsitution effect
  double pT;
  for(int i=0; i<X.n_cols; ++i){ //Matrix is column-major
    pT = p[i];
    dd += -d[i]*X.col(i)%(X.col(i)-2*pT-1)-2*d[i]*pT*pT;
  }
  X.each_row() -= 2*p;
  bv = geno*alpha;
  return Rcpp::List::create(Rcpp::Named("bv")=bv,
                            Rcpp::Named("dd")=dd,
                            Rcpp::Named("alpha")=alpha);
}

