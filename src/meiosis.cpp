// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
arma::Row<unsigned char> simGamHal2(const arma::Row<unsigned char>& chr1,
                                    const arma::Row<unsigned char>& chr2,
                                    const arma::vec& genMap){
  int nSites = chr1.n_elem;
  double genLen = genMap[nSites-1];
  arma::Row<unsigned char> gamete(nSites);
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

//' @title Diploid cross
//' @description Makes crosses between diploid individuals.
//' @param fGeno female genotypes
//' @param fPar female parents
//' @param mGeno male genotypes
//' @param mPar male parents
//' @param genMaps chromosome genetic maps
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
Rcpp::List cross2(const Rcpp::List& fGeno, arma::uvec fPar,
                  const Rcpp::List& mGeno, arma::uvec mPar,
                  const Rcpp::List& genMaps){
  // R to C++
  fPar = fPar-1;
  mPar = mPar-1;
  int nChr = fGeno.length();
  int nInd = fPar.n_elem;
  Rcpp::List geno(nChr);
  for(int i=0; i<nChr; ++i){
    arma::vec genMap = Rcpp::as<arma::vec>(genMaps[i]);
    Rcpp::List tmpOut(2);
    Rcpp::List tmpIn;
    arma::Mat<unsigned char> chr1;
    arma::Mat<unsigned char> chr2;
    //Make female gametes
    tmpIn = fGeno[i];
    chr1 = Rcpp::as<arma::Mat<unsigned char> >(tmpIn[0]);
    chr2 = Rcpp::as<arma::Mat<unsigned char> >(tmpIn[1]);
    int segSites = chr1.n_cols;
    arma::Mat<unsigned char> fGam(nInd,segSites);
    for(int j=0; j<nInd; ++j)
      fGam.row(j) = simGamHal2(chr1.row(fPar[j]),chr2.row(fPar[j]),genMap);
    tmpOut[0] = fGam;
    //Make male gametes
    tmpIn = mGeno[i];
    chr1 = Rcpp::as<arma::Mat<unsigned char> >(tmpIn[0]);
    chr2 = Rcpp::as<arma::Mat<unsigned char> >(tmpIn[1]);
    arma::Mat<unsigned char> mGam(nInd,segSites);
    for(int j=0; j<nInd; ++j)
      mGam.row(j) = simGamHal2(chr1.row(mPar[j]),chr2.row(mPar[j]),genMap);
    tmpOut[1] = mGam;
    //Return chromosomes
    geno[i] = tmpOut;
  }
  return geno;
}




