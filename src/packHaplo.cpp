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
    nInd = haplo.n_cols/ploidy;
  }
  arma::uword nLoci = haplo.n_cols;
  arma::Cube<unsigned char> output(nLoci,ploidy,nInd);
  arma::uword chr=0;
  for(arma::uword i=0; i<nHap; ++i){
    if(inbred){
      for(chr=0; chr<ploidy; ++chr){
        output.slice(i).col(chr) = haplo.row(i).t();
      }
    }else{
      output.slice(i/ploidy).col(chr) = haplo.row(i).t();
      ++chr;
      chr = chr%ploidy; 
    }
  }
  return output;
}

