#include "alphasimr.h"

// [[Rcpp::export]]
arma::Cube<unsigned char> packHaplo(arma::Mat<unsigned char>& haplo,
                                    arma::uword ploidy, bool inbred){
  arma::uword nHap = haplo.n_rows;
  arma::uword nInd;
  if(inbred){
    nInd = haplo.n_rows;
  }else{
    if(haplo.n_rows%ploidy != 0){
      Rcpp::stop("Number of rows not a factor of ploidy");
    }
    nInd = haplo.n_rows/ploidy;
  }
  arma::uword nLoci = haplo.n_cols;
  arma::uword nBins = nLoci/8;
  if((nLoci%8) > 0){
    ++nBins;
  }
  arma::Cube<unsigned char> output(nBins,ploidy,nInd);
  arma::uword chr=0;
  arma::uword locus;
  std::bitset<8> workBits;
  for(arma::uword i=0; i<nHap; ++i){
    if(inbred){
      for(chr=0; chr<ploidy; ++chr){
        locus = 0;
        
        // Fill in bins known to be complete
        if(nBins > 1){
          for(arma::uword j=0; j<(nBins-1); ++j){
            for(arma::uword k=0; k<8; ++k){
              workBits[k] = haplo(i,locus);
              ++locus;
            }
            output.slice(i).col(chr).row(j) = 
              toByte(workBits);
          }
        }
        
        // Fill in potentially incomplete bins
        for(arma::uword k=0; k<8; ++k){
          if(locus<nLoci){
            workBits[k] = haplo(i,locus);
            ++locus;
          }else{
            workBits[k] = 0;
          }
        }
        output.slice(i).col(chr).row(nBins-1) = 
          toByte(workBits);
      }
    }else{
      locus = 0;
      
      // Fill in bins known to be complete
      if(nBins > 1){
        for(arma::uword j=0; j<(nBins-1); ++j){
          for(arma::uword k=0; k<8; ++k){
            workBits[k] = haplo(i,locus);
            ++locus;
          }
          output.slice(i/ploidy).col(chr).row(j) = 
            toByte(workBits);
        }
      }
      
      // Fill in potentially incomplete bins
      for(arma::uword k=0; k<8; ++k){
        if(locus<nLoci){
          workBits[k] = haplo(i,locus);
          ++locus;
        }else{
          workBits[k] = 0;
        }
      }
      output.slice(i/ploidy).col(chr).row(nBins-1) = 
        toByte(workBits);
      ++chr;
      chr = chr%ploidy; 
    }
  }
  return output;
}
