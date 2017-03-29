// A stop-gap solution until MaCS is integrated

// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Reads haplotypes from AlphaFormatter as a cube of unsigned char
// [[Rcpp::export]]
arma::Cube<unsigned char> readAF(int nInd, int segSites, int ploidy,
                                 arma::uvec keep, bool inbred){
  int nLoci = keep.n_elem;
  // Calculated expected number of haplotypes
  int nHap;
  if(inbred){
    nHap = nInd;
  }else{
    nHap = nInd*ploidy;
  }
  arma::Cube<unsigned char> output(nLoci,ploidy,nInd);
  arma::Col<unsigned char> tmp(nLoci);
  // Read output from AlphaFormatter
  char sep = ' ';
  std::string iFileName = "MacsHaplotypes.txt";
  std::ifstream iFile(iFileName.c_str());
  std::string line;
  // Read rows
  int chr=0;
  for(int i=0; i<nHap; ++i){
    std::getline(iFile,line);
    std::stringstream lineStream(line);
    std::string cell;
    // Skip blank column
    std::getline(lineStream,cell,sep);
    // Read columns, saving data in uvec keep to tmp
    int k=0;
    int j=0;
    while(k<nLoci){
      std::getline(lineStream,cell,sep);
      if(j==keep(k)){
        tmp(k) = std::atof(cell.c_str());
        ++k;
      }
      ++j;
    }
    // Transfer tmp to output
    if(inbred){
      for(chr=0; chr<ploidy; ++chr){
        output.slice(i).col(chr) = tmp;
      }
    }else{
     output.slice(i/ploidy).col(chr) = tmp;
     ++chr;
      chr = chr%ploidy; 
    }
  }
  iFile.close();
  return output;
}
