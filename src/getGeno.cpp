// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const Rcpp::S4& pop, const Rcpp::S4& lociMap){
  int nInd = pop.slot("nInd");
  int nChr = pop.slot("nChr");
  int ploidy = pop.slot("ploidy");
  Rcpp::List geno = pop.slot("geno");
  arma::ivec lociPerChr = lociMap.slot("lociPerChr");
  arma::uvec lociLoc = lociMap.slot("lociLoc");
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr), arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(int i=0; i<nChr; ++i){
    loc1 = loc2+1;
    loc2 += lociPerChr[i];
    Rcpp::List chrGeno = geno[i];
    arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2))-1; //R to C++
    for(int j=0; j<ploidy; ++j)
      output.cols(loc1,loc2) +=
        Rcpp::as<arma::Mat<unsigned char> >(chrGeno[j]).cols(chrLociLoc);
  }
  return output;
}

//' @title Get SNP chip genotype data
//' 
//' @description Retrieves gentoype information from a population for a SNP chip.
//' 
//' @param pop an object of superclass 'Pop'
//' @param chip which chip
//' @param simParam an object of class 'SimParam'
//' 
//' @export
// [[Rcpp::export]]
arma::imat pullSnpGeno(const Rcpp::S4& pop, int chip,
                       const Rcpp::S4& simParam){
  Rcpp::List snpChips = simParam.slot("snpChips");
  Rcpp::S4 lociMap = snpChips[(chip-1)];
  return arma::conv_to<arma::imat>::from(getGeno(pop,lociMap));
}

//' @title Get QTL genotype data
//' 
//' @description Retrieves gentoype information from a population for a SNP chip.
//' 
//' @param pop an object of superclass 'Pop'
//' @param trait which trait
//' @param simParam an object of class 'SimParam'
//' 
//' @export
// [[Rcpp::export]]
arma::imat pullQtlGeno(const Rcpp::S4& pop, int trait,
                       const Rcpp::S4& simParam){
  Rcpp::List traits = simParam.slot("traits");
  Rcpp::S4 lociMap = traits[(trait-1)];
  return arma::conv_to<arma::imat>::from(getGeno(pop,lociMap));
}

// Converts geno matrix to dominance indicator matrix when ploidy=2
// Using imat instead of Mat<unsigned char> for increased speed
// [[Rcpp::export]]
arma::imat getDomGeno(const arma::Mat<unsigned char>& geno){
  arma::imat output(geno.n_rows,geno.n_cols);
  output = arma::conv_to<arma::imat>::from(geno);
  output -= 1;
  output = 1-arma::abs(output);
  return output;
}

// Calculates minor allele frequency on a single chromsomes for diploids
// [[Rcpp::export]]
arma::rowvec calcQ2(Rcpp::List& geno){
  arma::Mat<unsigned char> tmp;
  tmp = Rcpp::as<arma::Mat<unsigned char> >(geno[0]);
  tmp += Rcpp::as<arma::Mat<unsigned char> >(geno[1]);
  arma::mat X = arma::conv_to<arma::mat>::from(tmp);
  arma::rowvec output = arma::mean(X,0)/2;
  return 0.5-arma::abs(output-0.5);
}

