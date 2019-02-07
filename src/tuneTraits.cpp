// Functions for adjusting traits to desired initial variance
// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

//Tunes TraitA for desired varG
// [[Rcpp::export]]
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno,
                      arma::vec& addEff,
                      double varG, 
                      int ploidy, 
                      int nThreads){
  arma::vec gv(geno.n_cols, arma::fill::zeros);
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      gv(i) += xa(geno(j,i))*addEff(j);
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
                      bool useVarA,
                      int ploidy, 
                      int nThreads){
  double scale, obsVarG, obsVarA, intercept; //Output
  arma::uword nLoci = geno.n_rows;
  arma::uword nInd = geno.n_cols;
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  arma::vec genoMu(nLoci), alpha(nLoci);
  arma::vec gv(nInd,arma::fill::zeros), bv(nInd,arma::fill::zeros);
  arma::Mat<unsigned char> genoT = geno.t();
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<genoT.n_cols; ++i){
    arma::vec freq(ploidy+1,arma::fill::zeros);
    arma::vec gvLoc(ploidy+1);
    double gvMu;
    for(arma::uword j=0; j<genoT.n_rows; ++j){
      freq(genoT(j,i)) += 1;
    }
    freq = freq/accu(freq);
    genoMu(i) = accu(freq%x);
    gvLoc = xa*addEff(i)+xd*domEff(i);
    gvMu = accu(freq%gvLoc);
    
    alpha(i) = accu(freq%(gvLoc-gvMu)%(x-genoMu(i)))/
      accu(freq%(x-genoMu(i))%(x-genoMu(i)));
  }

  // Account for divide by zero
  alpha.elem(find_nonfinite(alpha)).zeros();
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      gv(i) += xa(geno(j,i))*addEff(j)+xd(geno(j,i))*domEff(j);
      bv(i) += (double(geno(j,i))-genoMu(j))*alpha(j);
    }
  }
  
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
