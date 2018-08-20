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

// Converts geno matrix to dominance indicator matrix when ploidy=2
// [[Rcpp::export]]
arma::Mat<unsigned char> getDomGeno(arma::Mat<unsigned char> geno){
  return geno.replace(2,0);
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

// Returns IBD haplotype data in a matrix of nInd*2 by nLoci
// [[Rcpp::export]]
Rcpp::IntegerMatrix getIbdHaplo(const Rcpp::IntegerMatrix & pedigree,
                                const Rcpp::List          & recHist,
                                const Rcpp::IntegerVector & lociPerChr) {
  int nInd = pedigree.nrow();
  int nChr = lociPerChr.size();
  int nLoc = sum(lociPerChr);
  int pId;
  Rcpp::IntegerMatrix output(nInd * 2, nLoc);
  Rcpp::IntegerVector gametes(nLoc); // all gametes passed from a parent to an individual
  for (int ind = 0; ind < nInd; ++ind) {
    // std::cout << ind + 1 << " " << pedigree(ind, 0) << " " << pedigree(ind, 1) << "\n";
    Rcpp::List recHistInd = recHist[ind];
    for (int par = 0; par < 2; ++par) {
      pId = pedigree(ind, par); // note AlphaSimR has mother/father as the first/second parent
      if (pId == 0) {
        gametes = Rcpp::rep(2 * ind + par + 1, nLoc); // first gamete is maternal, second is paternal
      } else {
        int chrOrigin = 0;
        for (int chr = 0; chr < nChr; ++chr) {
          if (lociPerChr[chr] > 0) {
            Rcpp::List recHistIndChr = recHistInd[chr];
            Rcpp::IntegerMatrix recHistIndChrPar = recHistIndChr[par];
            int nSeg = recHistIndChrPar.nrow();
            // std::cout << chr + 1 << " " << par + 1 << " " << nSeg << "\n";
            for (int seg = 0; seg < nSeg; ++seg) {
              int source = recHistIndChrPar(seg, 0);
              int start = chrOrigin + recHistIndChrPar(seg, 1);
              int stop;
              if (seg < (nSeg - 1)) {
                stop = chrOrigin + recHistIndChrPar(seg + 1, 1) - 1;
              } else {
                stop = chrOrigin + lociPerChr[chr];
              }
              // std::cout << start << " " << stop << "\n";
              for (int loc = start - 1; loc < stop; ++loc) {
                // note pId    is 1:n, hence (pid    - 1)
                // note source is 1:2, hence (source - 1)
                gametes[loc] = output(2 * (pId - 1) + (source - 1), loc);
              }
            }
          }
          chrOrigin = chrOrigin + lociPerChr[chr];
        }
      }
      output(2 * ind + par, Rcpp::_) = gametes;
    }
  }
  return output;
}

// Returns IBD haplotype data in a matrix of nInd*2 by nLoci
// [[Rcpp::export]]
Rcpp::IntegerMatrix getIbdHaplo2(const Rcpp::List          & ibdRecHist,
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
