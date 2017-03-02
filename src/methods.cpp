// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @title Get genotype data
//' @description Retrieves allele dossages for requested loci.
//' @param geno raw genotype data
//' @param nInd number of individuals in population
//' @param nChr number of chromosomes
//' @param ploidy ploidy of species
//' @param lociPerChr number of loci per chromosome
//' @param lociLoc physical position of loci
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(Rcpp::List& geno, int nInd, int nChr, int ploidy,
                                 arma::ivec& lociPerChr, arma::uvec& lociLoc){
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
