// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @title Get genotype
//' 
//' @description Retrieves gentoype information from a population for a trait or SNP chip.
//' 
//' @param pop an object of superclass 'Pop'
//' @param lociMap an object of superclass 'LociMap'
//' 
//' @export
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



