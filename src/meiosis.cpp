// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"
#include <random>

// Simulates crossing over during meiosis
// May be extended later to include higher ploidy levels

// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns -1 if value is smaller than the values of x
// Returns last element if value is greater than values of x
// Set left to the smallest value of the interval to search
int intervalSearch(arma::vec x, double value, int left=0){
  // Check if crossover is before beginning
  if(x[left]>value){
    // Return error
    return -1;
  }
  int end = x.n_elem-1;
  // Check if crossover is at or past end
  if(x[end]<=value){
    return end;
  }
  // Perform search
  int right = end;
  while((right-left)>1){ // Interval can be decreased
    int middle = (left + right) / 2;
    if (x[middle] == value){
      left = middle;
      // Check if at the end of the vector
      if(left<end){
        // Check for identical values to the right
        while(x[left+1]==value){
          left += 1;
          if(left==end){
            break;
          }
        }
      }
      break;
    } else if (x[middle] > value){
      right = middle;
    }else{
      left = middle;
    }
  }
  return left;
}

//Simulates a gamete using Haldane's model for crossovers, ploidy=2
arma::Col<unsigned char> bivalent(const arma::Col<unsigned char>& chr1,
                                  const arma::Col<unsigned char>& chr2,
                                  const arma::vec& genMap){
  int nSites = chr1.n_elem;
  double genLen = genMap(nSites-1);
  arma::Col<unsigned char> gamete(nSites);
  int nCO = Rcpp::rpois(1, genLen)(0);
  if(nCO==0){
    // No CO, randomly pick a chromosome
    if(Rcpp::rbinom(1,1,0.5)(0)){
      gamete = chr1;
    }else{
      gamete = chr2;
    }
  }else{
    // COs present, Create CO locations
    arma::vec posCO(nCO,arma::fill::randu);
    posCO = posCO*genLen;
    posCO = arma::sort(posCO);
    int startPos = 0;
    int endPos;
    // Randomly pick starting chromosome and fill first segSite
    int readChr = Rcpp::rbinom(1,1,0.5)(0);
    if(readChr){
      gamete(0) = chr1(0);
    }else{
      gamete(0) = chr2(0);
    }
    for(int i=0; i<nCO; ++i){
      endPos = intervalSearch(genMap,posCO[i],startPos);
      // Fill gamete
      if(endPos>startPos){ // Check for double crossovers
        // Fill in segSites
        if(readChr){
          gamete(arma::span(startPos+1,endPos)) = chr1(arma::span(startPos+1,endPos));
        }else{
          gamete(arma::span(startPos+1,endPos)) = chr2(arma::span(startPos+1,endPos));
        }
      }
      startPos = endPos;
      // Switch chromosome
      ++readChr;
      readChr = Rcpp::rbinom(1,1,0.5)(0);
    }
    // Fill in last segSites if needed
    if(endPos<(nSites-1)){
      if(readChr){
        gamete(arma::span(endPos+1,nSites-1)) = chr1(arma::span(endPos+1,nSites-1));
      }else{
        gamete(arma::span(endPos+1,nSites-1)) = chr2(arma::span(endPos+1,nSites-1));
      }
    }
  }
  return gamete;
}

// Makes crosses between diploid individuals.
// fGeno: female genotypes
// fPar: female parents
// mGeno: male genotypes
// mPar: male parents
// genMaps: chromosome genetic maps
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > cross2(
    const arma::field<arma::Cube<unsigned char> >& fGeno, 
    arma::uvec fPar,
    const arma::field<arma::Cube<unsigned char> >& mGeno, 
    arma::uvec mPar,
    const arma::field<arma::vec>& genMaps){
  fPar -= 1; // R to C++
  mPar -= 1; // R to C++
  int nChr = fGeno.n_elem;
  int nInd = fPar.n_elem;
  //Output data
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  //Loop through chromosomes
  for(int chr=0; chr<nChr; ++chr){
    int segSites = fGeno(chr).n_rows;
    arma::Cube<unsigned char> tmpGeno(segSites,2,nInd);
    //Loop through individuals
    for(int ind=0; ind<nInd; ++ind){
      //Female gamete
      tmpGeno.slice(ind).col(0) = 
        bivalent(fGeno(chr).slice(fPar(ind)).col(0),
                 fGeno(chr).slice(fPar(ind)).col(1),
                 genMaps(chr));
      //Male gamete
      tmpGeno.slice(ind).col(1) = 
        bivalent(mGeno(chr).slice(mPar(ind)).col(0),
                 mGeno(chr).slice(mPar(ind)).col(1),
                 genMaps(chr));
    } //End individual loop
    geno(chr) = tmpGeno;
  } //End chromosome loop
  return geno;
}

// Creates DH lines from diploid individuals
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > createDH2(
    const arma::field<arma::Cube<unsigned char> >& geno, 
    int nDH, const arma::field<arma::vec>& genMaps){
  int nChr = geno.n_elem;
  int nInd = geno(0).n_slices;
  //Output data
  arma::field<arma::Cube<unsigned char> > output(nChr);
  for(int chr=0; chr<nChr; ++chr){ //Chromosome loop
    int segSites = geno(chr).n_rows;
    arma::Cube<unsigned char> tmp(segSites,2,nInd*nDH);
    for(int ind=0; ind<nInd; ++ind){ //Individual loop
      for(int i=0; i<nDH; ++i){ //nDH loop
        arma::Col<unsigned char> gamete = 
          bivalent(geno(chr).slice(ind).col(0),
                   geno(chr).slice(ind).col(1),
                   genMaps(chr));
        for(int j=0; j<2; ++j){ //ploidy loop
          tmp.slice(i+ind*nDH).col(j) = gamete;
        } //End ploidy loop
      } //End nDH loop
    } //End individual loop
    output(chr) = tmp;
  } //End chromosome loop
  return output;
}

// Makes crosses between diploid individuals.
// fGeno: female genotypes
// fPar: female parents
// mGeno: male genotypes
// mPar: male parents
// genMaps: chromosome genetic maps
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > crossPedigree(
    const arma::field<arma::Cube<unsigned char> >& founders, 
    arma::uvec fPar,
    arma::uvec mPar,
    const arma::field<arma::vec>& genMaps){
  fPar -= 1; // R to C++
  mPar -= 1; // R to C++
  int nChr = founders.n_elem;
  int nInd = fPar.n_elem;
  
  typedef std::minstd_rand G;
  G g;
  typedef std::uniform_int_distribution<> D;
  D d(0,founders(0).n_slices-1);
  
  //Output data
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  //Loop through chromosomes
  for(int chr=0; chr<nChr; ++chr){
    int segSites = founders(chr).n_rows;
    arma::Cube<unsigned char> tmpGeno(segSites,2,nInd);
    
    //Loop through individuals
    for(int ind=0; ind<nInd; ++ind){
      if (fPar(ind) == -1){
        //Female gamete
        tmpGeno.slice(ind).col(0) = 
          bivalent(founders(chr).slice(d(g)).col(0),
                   founders(chr).slice(d(g)).col(1),
                   genMaps(chr));
      }
      else {
        //Female gamete
        tmpGeno.slice(ind).col(0) = 
          bivalent(tmpGeno.slice(fPar(ind)).col(0),
                   tmpGeno.slice(fPar(ind)).col(1),
                   genMaps(chr));           
      }
      if (mPar(ind) == -1) {
        //Male gamete
        tmpGeno.slice(ind).col(1) = 
          bivalent(founders(chr).slice(d(g)).col(0),
                   founders(chr).slice(d(g)).col(1),
                   genMaps(chr));
      }
      else
      {
        //Male gamete
        tmpGeno.slice(ind).col(1) = 
          bivalent(tmpGeno.slice(mPar(ind)).col(0),
                   tmpGeno.slice(mPar(ind)).col(1),
                   genMaps(chr));
      }
    } //End individual loop
    geno(chr) = tmpGeno;
  } //End chromosome loop
  return geno;
}
