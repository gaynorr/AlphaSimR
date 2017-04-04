// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

/*
 * Genotype data is stored in a field of cubes.
 * The field has length equal to nChr
 * Each cube has dimensions nLoci by ploidy by nInd
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::ivec& lociPerChr,
                                 arma::uvec lociLoc){
  // R to C++ index correction
  lociLoc -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr), arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(int i=0; i<nChr; ++i){
    // Get loci locations
    loc1 = loc2+1;
    loc2 += lociPerChr[i];
    arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
    // Get chromsome genotype
    arma::Mat<unsigned char> tmp;
    tmp = arma::sum(geno(i),1);
    // Assign genotypes to output matrix
    output.cols(loc1,loc2) = (tmp.rows(chrLociLoc)).t();
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

// Returns haplotype data in a matrix of nInd*ploidy by nLoci
// [[Rcpp::export]]
arma::Mat<unsigned char> getHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                  const arma::ivec& lociPerChr,
                                  arma::uvec lociLoc){
  // R to C++ index correction
  lociLoc -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  int ploidy = geno(0).n_cols;
  int nLoci = geno(0).n_rows;
  arma::Mat<unsigned char> output(nInd*ploidy,arma::sum(lociPerChr), arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  arma::Mat<unsigned char> tmp(nLoci,ploidy);
  // Get chromosome data
  for(int i=0; i<nChr; ++i){
    // Get loci locations
    loc1 = loc2+1;
    loc2 += lociPerChr[i];
    arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
    // Get individual data
    for(int ind=0; ind<nInd; ++ind){
      output(arma::span(ind*ploidy,(ind+1)*ploidy-1),
             arma::span(loc1,loc2)) = 
        (geno(i).slice(ind).rows(chrLociLoc)).t();
    }
  }
  return output;
}