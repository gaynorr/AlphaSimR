// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

/*
 * Genotype data is stored in a field of cubes.
 * The field has length equal to nChr
 * Each cube has dimensions nLoci by ploidy by nInd
 * Output return with dimensions nInd by nLoci
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::ivec& lociPerChr,
                                 arma::uvec lociLoc){
  // R to C++ index correction
  lociLoc -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr));
  int loc1;
  int loc2 = -1;
  for(int i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      // Get chromsome genotype
      arma::Mat<unsigned char> tmp;
      tmp = arma::sum(geno(i),1);
      // Assign genotypes to output matrix
      output.cols(loc1,loc2) = (tmp.rows(chrLociLoc)).t();
    }
  }
  return output;
}

/*
 * Genotype data is stored in a field of cubes.
 * The field has length equal to nChr
 * Each cube has dimensions nLoci by ploidy by nInd
 * Output return with dimensions nLoci by nInd
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGenoT(const arma::field<arma::Cube<unsigned char> >& geno, 
                                  const arma::ivec& lociPerChr,
                                  arma::uvec lociLoc){
  // R to C++ index correction
  lociLoc -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  arma::Mat<unsigned char> output(arma::sum(lociPerChr),nInd);
  int loc1;
  int loc2 = -1;
  for(int i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      // Get chromsome genotype
      arma::Mat<unsigned char> tmp;
      tmp = arma::sum(geno(i),1);
      // Assign genotypes to output matrix
      output.rows(loc1,loc2) = tmp.rows(chrLociLoc);
    }
  }
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
  arma::Mat<unsigned char> output(nInd*ploidy,arma::sum(lociPerChr));
  int loc1;
  int loc2 = -1;
  // Get chromosome data
  for(int i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      // Get individual data
      for(int ind=0; ind<nInd; ++ind){
        output(arma::span(ind*ploidy,(ind+1)*ploidy-1),
               arma::span(loc1,loc2)) = 
                 (geno(i).slice(ind).rows(chrLociLoc)).t();
      }
    }
  }
  return output;
}

// Returns haplotype data in a matrix of nInd by nLoci for a single
// chromosome group. i.e. just female or male chromosomes for diploids
// [[Rcpp::export]]
arma::Mat<unsigned char> getOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                     const arma::ivec& lociPerChr,
                                     arma::uvec lociLoc, int haplo){
  // R to C++ index correction
  lociLoc -= 1;
  haplo -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr));
  int loc1;
  int loc2 = -1;
  arma::uvec colSel(1);
  colSel(0) = haplo;
  // Get chromosome data
  for(int i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      // Get individual data
      for(int ind=0; ind<nInd; ++ind){
        output(ind,arma::span(loc1,loc2)) = 
          (geno(i).slice(ind).submat(chrLociLoc,colSel)).t();
      }
    }
  }
  return output;
}

// Returns haplotype data in a matrix of nLoci by nInd for a single
// chromosome group. i.e. just female or male chromosomes for diploids
// [[Rcpp::export]]
arma::Mat<unsigned char> getOneHaploT(const arma::field<arma::Cube<unsigned char> >& geno, 
                                      const arma::ivec& lociPerChr,
                                      arma::uvec lociLoc, int haplo){
  // R to C++ index correction
  lociLoc -= 1;
  haplo -= 1;
  
  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  arma::Mat<unsigned char> output(arma::sum(lociPerChr),nInd);
  int loc1;
  int loc2 = -1;
  arma::uvec colSel(1);
  colSel(0) = haplo;
  // Get chromosome data
  for(int i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      // Get individual data
      for(int ind=0; ind<nInd; ++ind){
        output(arma::span(loc1,loc2),ind) = 
          geno(i).slice(ind).submat(chrLociLoc,colSel);
      }
    }
  }
  return output;
}

// Returns IBD haplotype data in a matrix of nInd*2 by nLoci
// [[Rcpp::export]]
Rcpp::IntegerMatrix getIbdHaplo(const Rcpp::List          & ibdRecHist,
                                const Rcpp::IntegerVector & individuals,
                                const Rcpp::IntegerVector & nLociPerChr) {
  int nIndSet = individuals.size();
  int nChr = nLociPerChr.size();
  int nLoc = sum(nLociPerChr);
  Rcpp::IntegerMatrix output(nIndSet * 2, nLoc);
  for (int indSet = 0; indSet < nIndSet; ++indSet) {
    // std::cout << indSet + 1 << "\n";
    Rcpp::List ibdRecHistInd = ibdRecHist(individuals(indSet) - 1);
    for (int par = 0; par < 2; ++par) {
      int chrOrigin = 0;
      for (int chr = 0; chr < nChr; ++chr) {
        if (nLociPerChr(chr) > 0) {
          Rcpp::List ibdRecHistIndChr = ibdRecHistInd(chr);
          Rcpp::IntegerMatrix ibdRecHistIndChrPar = ibdRecHistIndChr(par);
          int nSeg = ibdRecHistIndChrPar.nrow();
          // std::cout << chr + 1 << " " << par + 1 << " " << nSeg << "\n";
          for (int seg = 0; seg < nSeg; ++seg) {
            int source = ibdRecHistIndChrPar(seg, 0);
            int start = chrOrigin + ibdRecHistIndChrPar(seg, 1);
            int stop;
            if (seg < (nSeg - 1)) {
              stop = chrOrigin + ibdRecHistIndChrPar(seg + 1, 1) - 1;
            } else {
              stop = chrOrigin + nLociPerChr[chr];
            }
            // std::cout << start << " " << stop << "\n";
            for (int loc = start - 1; loc < stop; ++loc) {
              output(2 * indSet + par, loc) = source;
            }
          }
        }
        chrOrigin = chrOrigin + nLociPerChr(chr);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
void writeGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
               const arma::ivec& lociPerChr,
               arma::uvec lociLoc,
               Rcpp::String filePath){
  arma::Mat<unsigned char> output;
  output = getGeno(geno,lociPerChr,lociLoc);
  std::ofstream outFile;
  outFile.open(filePath, std::ios_base::app);
  output.save(outFile,arma::raw_ascii);
  outFile.close();
}

// [[Rcpp::export]]
void writeOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                   const arma::ivec& lociPerChr, 
                   arma::uvec lociLoc, int haplo,
                   Rcpp::String filePath){
  arma::Mat<unsigned char> output;
  output = getOneHaplo(geno,lociPerChr,lociLoc,haplo);
  std::ofstream outFile;
  outFile.open(filePath, std::ios_base::app);
  output.save(outFile,arma::raw_ascii);
  outFile.close();
}
