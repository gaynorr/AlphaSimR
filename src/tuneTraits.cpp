// Functions for adjusting traits to desired initial variance
// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

//Tunes TraitA for desired varG
// [[Rcpp::export]]
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      double varG, int nThreads){
  arma::vec gv  = calcGvA(geno, addEff, 0, nThreads);
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
                      bool useVarA, 
                      int nThreads){
  int nLoci = geno.n_rows;
  int nInd = geno.n_cols;
  arma::vec p(nLoci), q(nLoci), alpha(nLoci);
  arma::vec gv(nInd,arma::fill::zeros), bv(nInd,arma::fill::zeros);
  arma::Mat<unsigned char> genoT = geno.t();
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(int i=0; i<nLoci; ++i){
    arma::vec genoFreq(3,arma::fill::zeros);
    for(int j=0; j<nInd; ++j){
      genoFreq(genoT(j,i)) += 1;
    }
    p(i) = (genoFreq(2)+0.5*genoFreq(1))/accu(genoFreq);
    q(i) = 1-p(i);
    // 1-observed(het)/expect(het)
    if((p(i)>0.999999999) | (p(i)<0.000000001)){
      // Locus is fixed, no viable regression
      alpha(i) = 0;
    }else{
      double F = 1-(genoFreq(1)/accu(genoFreq))/(2*p(i)*q(i));
      if(F<-0.999999999){
        // Only heterozygotes, no viable regression
        alpha(i) = 0;
      }else{
        // a+d(q-p)(1-F)/(1+F)
        alpha(i) = addEff(i)+domEff(i)*(q(i)-p(i))*(1-F)/(1+F);
      }
    }
  }
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      gv(i) += geno(j,i)*addEff(j)+(1-abs(int(geno(j,i))-1))*domEff(j);
      bv(i) += (double(geno(j,i))-2*p(j))*alpha(j);
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
