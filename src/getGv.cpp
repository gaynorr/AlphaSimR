// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates genetic values for a trait with only additive effects
arma::vec calcGvA(const arma::Mat<unsigned char>& geno,
                  const arma::vec& a, double intercept){
  arma::vec output = geno*a;
  return output+intercept;
}

// Calculates genetic values for a trait with additive and dominance effects
arma::vec calcGvAD(const arma::Mat<unsigned char>& geno,
                   const arma::vec& a, const arma::vec& d,
                   double intercept){
  arma::vec output = geno*a+getDomGeno(geno)*d;
  return output+intercept;
}

// Retrieves genetic values for TraitA
// A wrapper for accessing calcGvA
// [[Rcpp::export]]
arma::vec getGvA(const Rcpp::S4& trait, const Rcpp::S4& pop){
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  return calcGvA(geno, a, intercept);
}

// Retrieves genetic value for TraitAD
// A wrapper for accessing calcGvAD
// [[Rcpp::export]]
arma::vec getGvAD(const Rcpp::S4& trait, const Rcpp::S4& pop){
  arma::vec a = trait.slot("addEff");
  arma::vec d = trait.slot("domEff");
  double intercept = trait.slot("intercept");
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  return calcGvAD(geno, a, d, intercept);
}

// A calculates breeding values, dominance deviations and allele
// subsitution effects. Used for calculating additive and dominance
// genetic variances. Only works for ploidy=2.
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait, const Rcpp::S4& pop){
  arma::vec a = trait.slot("addEff");
  arma::vec d(a.n_elem);
  if(trait.hasSlot("domEff")){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }else{
    d.zeros();
  }
  int nInd = pop.slot("nInd");
  arma::vec bv(nInd);
  arma::vec dd(nInd,arma::fill::zeros);
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  arma::mat X = arma::conv_to<arma::mat>::from(geno);
  arma::rowvec p = arma::mean(X,0)/2.0;
  arma::vec alpha = a+d%(1-2*p.t()); //allele subsitution effect
  double pT;
  for(int i=0; i<X.n_cols; ++i){ //Matrix is column-major
    pT = p(i);
    dd += -d(i)*X.col(i)%(X.col(i)-2*pT-1)-2*d(i)*pT*pT;
  }
  X.each_row() -= 2*p;
  bv = X*alpha;
  arma::vec q = 1-p.t();
  double genicVarA = 2.0*arma::sum(p.t()%q%arma::square(alpha));
  double genicVarD = 4.0*arma::sum(arma::square(p.t())%arma::square(q)%arma::square(d));
  return Rcpp::List::create(Rcpp::Named("bv")=bv,
                            Rcpp::Named("dd")=dd,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("genicVarA")=genicVarA,
                            Rcpp::Named("genicVarD")=genicVarD);
}

