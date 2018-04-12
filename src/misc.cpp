// These function are called by R, but are not listed in the package namespace
// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates population variance
//' @title Population variance
//' 
//' @description
//' Calculates the population variance matrix as 
//' opposed to the sample variance matrix calculated 
//' by \code{\link{var}}. i.e. divides by n instead 
//' of n-1
//' 
//' @param X an n by m matrix
//' 
//' @return an m by m variance-covariance matrix
//' 
//' @export
// [[Rcpp::export]]
arma::mat popVar(const arma::mat& X) {
  return arma::cov(X,1);
}

// Merges geno objects, i.e. fields containing cubes of unsigned char
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > mergeGeno(
    const arma::field<arma::Cube<unsigned char> >& x, 
    const arma::field<arma::Cube<unsigned char> >& y){
  int nChr = x.n_elem;
  arma::field<arma::Cube<unsigned char> > z(nChr);
  for(arma::uword i=0; i<nChr; ++i){
    z(i) = arma::join_slices(x(i),y(i));
  }
  return z;
}

// Merges multiple geno objects contained a list of Class-Pop
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > mergeMultGeno(Rcpp::List& popList,
                                                      arma::uvec nInd,
                                                      arma::uvec nLoci,
                                                      arma::uword ploidy){
  arma::field<arma::Cube<unsigned char> > output(nLoci.n_elem);
  arma::uword nTot = sum(nInd);
  arma::uword nPop = nInd.n_elem;
  // Allocate output
  for(arma::uword chr=0; chr<nLoci.n_elem; ++chr){
    output(chr).set_size(nLoci(chr),ploidy,nTot);
  }
  // Add individual genotypes
  arma::uword startInd=0, endInd=0;
  for(arma::uword i=0; i<nPop; ++i){
    if(nInd(i)>0){
      endInd += nInd(i)-1;
      Rcpp::S4 pop = popList[i];
      arma::field<arma::Cube<unsigned char> >geno = pop.slot("geno");
      for(arma::uword chr=0; chr<nLoci.n_elem; ++chr){
        output(chr).slices(startInd,endInd) = geno(chr);
      }
      startInd += nInd(i);
      endInd = startInd;
    }
  }
  return output;
}

// Merges a list of integer matrices
// [[Rcpp::export]]
arma::Mat<int> mergeMultIntMat(const arma::field<arma::Mat<int> >& X,
                               arma::uvec nRow,
                               arma::uword nCol){
  arma::Mat<int> output(sum(nRow),nCol);
  arma::uword start=0, end=0;
  for(arma::uword i=0; i<nRow.n_elem; i++){
    if(nRow(i)>0){
      end += nRow(i)-1;
    }
    output.rows(start,end) = X(i);
    start += nRow(i);
    end = start;
  }
  return output;
}

// Calculates allele frequency on a single chromsome
// Requires bi-allelic markers, but works for any ploidy
// [[Rcpp::export]]
arma::vec calcChrFreq(const arma::Cube<unsigned char>& geno){
  int ploidy = geno.n_cols;
  arma::Mat<unsigned char> tmp = arma::sum(geno,1);
  arma::vec output = arma::mean(arma::conv_to<arma::mat>::from(tmp),
                                1)/ploidy;
  return output;
}

// [[Rcpp::export]]
arma::Mat<int> convToImat(const arma::Mat<unsigned char>& X){
  return arma::conv_to<arma::Mat<int> >::from(X);
}

// Linear index functions for upper triangle of a square
// matrix without the diagonal
// From: https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix

// Find mapping index
// i = row of matrix
// j = column of matrix
// n = dimension of matrix (i.e. row/column length)
arma::uword mapIndex(arma::uword i, arma::uword j,
                     arma::uword n){
  return (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j-i-1;
}

// Find row given mapping index
// k = mapping index
// n = dimension of matrix
arma::uword mapRow(arma::uword k, arma::uword n){
  return n-2-static_cast<arma::uword>(sqrt(-8*double(k) + 4*double(n)*(double(n)-1)-7)/2-0.5);
}

// Find column given mapping index
// k = mapping index
// n = dimension of matrix
arma::uword mapCol(arma::uword k, arma::uword n){
  arma::uword i;
  i = mapRow(k,n);
  return k+i+1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
}


// Randomly samples integers without replacement
// n number of integers to return
// N number of integers to sample from
// Returns an integer vector of length n with values ranging from 0 to N-1
// From: https://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
// Reportedly from: Algorithm 3.4.2S of Knuth's book Seminumeric Algorithms
arma::Col<arma::uword> sampleInt(arma::uword n, arma::uword N){
  arma::uword t = 0;
  arma::uword m = 0;
  Rcpp::NumericVector u(1);
  arma::Col<arma::uword> samples(n);
  while(m<n){
    u = Rcpp::runif(1);
    if(double(N-t)*u(0) >= double(n-m)){
      ++t;
    }else{
      samples(m) = t;
      ++t;
      ++m;
    }
  }
  return samples;
}

// Samples random pairs without replacement from all possible combinations
// nLevel1 = number of levels for the first column
// nLevel2 = number of levels for the second column
// n = number of combinations to sample
// If n is larger than the total number of possible combinations (N), 
// then only n%N combinations are sampled and the rest are systematically assigned
// Returns an integer matrix with the sampled levels for each column
// Values in column 1 range from 1 to nLevel1
// Values in column 2 range from 1 to nLevel2
// [[Rcpp::export]]
arma::Mat<arma::uword> sampAllComb(arma::uword nLevel1, arma::uword nLevel2, 
                                   arma::uword n){
  arma::uword N = nLevel1*nLevel2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::Col<arma::uword> samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::Mat<arma::uword> output(n,2);
  for(arma::uword  i=0; i<n; ++i){
    output(i,0) = samples(i)/nLevel2;
    output(i,1) = samples(i)%nLevel2;
  }
  if(fullComb>0){
    arma::Mat<arma::uword> tmp(N*fullComb,2);
    arma::uword i;
    for(arma::uword j=0; j<(N*fullComb); ++j){
      i = j%N;
      tmp(j,0) = i/nLevel2;
      tmp(j,1) = i%nLevel2;
    }
    output = arma::join_cols(output,tmp);
  }
  // C++ to R
  output += 1;
  return output;
}

// Samples random pairs without replacement from all half-diallel combinations
// nLevel = number of levels (number of individuals)
// n = number of combinations to sample
// If n is larger than the total number of possible combinations (N), 
// then only n%N combinations are sampled and the rest are systematically assigned
// Returns an integer matrix with the sampled levels for each combination
// Returned values range from 1 to nLevel
// [[Rcpp::export]]
arma::Mat<arma::uword> sampHalfDialComb(arma::uword nLevel, arma::uword n){
  arma::uword N = nLevel*(nLevel-1)/2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::Col<arma::uword> samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::Mat<arma::uword> output(n,2);
  for(arma::uword i=0; i<n; ++i){
    output(i,0) = mapRow(samples(i),nLevel);
    output(i,1) = mapCol(samples(i),nLevel);
  }
  if(fullComb>0){
    arma::Mat<arma::uword> tmp(N*fullComb,2);
    arma::uword i;
    for(arma::uword j=0; j<(N*fullComb); ++j){
      i = j%N;
      tmp(j,0) = mapRow(i,nLevel);
      tmp(j,1) = mapCol(i,nLevel);
    }
    output = arma::join_cols(output,tmp);
  }
  // C++ to R
  output += 1;
  return output;
}

// [[Rcpp::export]]
arma::mat calcCoef(arma::mat& X, arma::mat& Y){
  return arma::solve(X,Y);
}
