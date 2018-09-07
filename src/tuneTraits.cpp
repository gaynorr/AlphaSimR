// Functions for adjusting traits to desired initial variance
// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

//Tunes TraitA for desired varG
// [[Rcpp::export]]
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      double varG){
  arma::vec gv(geno.n_rows,arma::fill::zeros);
  for(arma::uword j=0; j<geno.n_cols; ++j){
    for(arma::uword i=0; i<geno.n_rows; ++i){
      gv(i) += geno(i,j)*addEff(j);
    }
  }
  double scale = sqrt(varG)/stddev(gv,1);
  double intercept = mean(gv*scale);
  return Rcpp::List::create(Rcpp::Named("scale")=scale,
                            Rcpp::Named("intercept")=intercept);
}

//Tunes TraitAD for desired varG (or varA), tuneTraitA used for inbreds
// [[Rcpp::export]]
Rcpp::List tuneTraitAD(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      arma::vec& domEff,
                      double varG, 
                      bool useVarA){
    
  arma::vec p(geno.n_cols);
  for(arma::uword i=0; i<geno.n_cols; ++i){
    p(i) = mean(arma::conv_to<arma::vec>::from(geno.col(i)))/2;
  }
  arma::vec alpha = addEff+domEff%(1-2*p);
  arma::vec gv(geno.n_rows,arma::fill::zeros), bv(geno.n_rows,arma::fill::zeros);
  for(arma::uword j=0; j<geno.n_cols; ++j){
    for(arma::uword i=0; i<geno.n_rows; ++i){
      gv(i) += geno(i,j)*addEff(j)+(1-abs(int(geno(i,j))-1))*domEff(j);
      bv(i) += geno(i,j)*alpha(j);
    }
  }
  double scale, obsVarG, obsVarA, intercept;
  if(useVarA){
    scale = sqrt(varG)/stddev(bv,1);
    obsVarA = varG;
    gv *= scale;
    intercept = mean(gv);
    obsVarG = arma::var(gv,1);
  }else{
    scale = sqrt(varG)/stddev(gv,1);
    obsVarG = varG;
    intercept = mean(gv*scale);
    obsVarA = arma::var(bv*scale,1);
  }
  return Rcpp::List::create(Rcpp::Named("scale")=scale,
                            Rcpp::Named("intercept")=intercept,
                            Rcpp::Named("varA")=obsVarA,
                            Rcpp::Named("varG")=obsVarG);
}
