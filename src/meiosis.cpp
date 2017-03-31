// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

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
  int nCO = Rcpp::rpois(1, genLen)[0];
  if(nCO==0){
    // No CO, randomly pick a chromosome
    if(std::rand()%2){
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
    int readChr = std::rand()%2;
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
      readChr = readChr%2;
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




