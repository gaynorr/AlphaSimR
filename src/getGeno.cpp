// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
 * Retrieves allele dossages from requested loci
 * geno: raw genotype data
 * nInd: number of individuals in population
 * nChr: number of chromosomes
 * ploidy: ploidy of species
 * lociPerChr: number of loci per chromosome
 * lociLoc: physical position of loci
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const Rcpp::S4& pop, const arma::ivec& lociPerChr, 
                                 const arma::uvec& lociLoc){
  int nInd = pop.slot("nInd");
  int nChr = pop.slot("nChr");
  int ploidy = pop.slot("ploidy");
  Rcpp::List geno = pop.slot("geno");
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



