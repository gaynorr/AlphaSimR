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
  for(int i=0; i<nChr; ++i){
    z(i) = arma::join_slices(x(i),y(i));
  }
  return z;
}

// Calculates minor allele frequency on a single chromsome
// Requires bi-allelic markers, but works for any ploidy
// [[Rcpp::export]]
arma::vec calcChrMinorFreq(const arma::Cube<unsigned char>& geno,
                           int ploidy){
  arma::Mat<unsigned char> tmp = arma::sum(geno,1);
  arma::vec output = arma::mean(arma::conv_to<arma::mat>::from(tmp),
                                1)/ploidy;
  return 0.5-arma::abs(output-0.5);
}

// [[Rcpp::export]]
arma::imat convToImat(const arma::Mat<unsigned char>& X){
  return arma::conv_to<arma::imat>::from(X);
}

// Linear index functions for upper triangle of a square
// matrix without the diagonal
// From: https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix

// Find mapping index
// i = row of matrix
// j = column of matrix
// n = dimension of matrix (i.e. row/column length)
long long int mapIndex(long long int i, long long int j,
                       long long int n){
  return (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j-i-1;
}

// Find row given mapping index
// k = mapping index
// n = dimension of matrix
long long int mapRow(long long int k, long long int n){
  return n-2-static_cast<long long int>(sqrt(-8*k + 4*n*(n-1)-7)/2-0.5);
}

// Find column given mapping index
// k = mapping index
// n = dimension of matrix
long long int mapCol(long long int k, long long int n){
  long long int i;
  i = mapRow(k,n);
  return k+i+1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
}


// Randomly samples integers without replacement
// n number of integers to return
// N number of integers to sample from
// Returns an integer vector of length n with values ranging from 0 to N-1
// From: https://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
// Reportedly from: Algorithm 3.4.2S of Knuth's book Seminumeric Algorithms
arma::ivec sampleInt(long long int n, long long int N){
  long long int t = 0;
  long long int m = 0;
  Rcpp::NumericVector u(1);
  arma::ivec samples(n);
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
arma::imat sampAllComb(long long int nLevel1, long long int nLevel2, 
                         long long int n){
  long long int N = nLevel1*nLevel2;
  long long int fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::ivec samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::imat output(n,2);
  for(long long int i=0; i<n; ++i){
    output(i,0) = samples(i)/nLevel2;
    output(i,1) = samples(i)%nLevel2;
  }
  if(fullComb>0){
    arma::imat tmp(N*fullComb,2);
    long long int i;
    for(long long int j=0; j<(N*fullComb); ++j){
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
arma::imat sampHalfDialComb(long long int nLevel, long long int n){
  long long int N = nLevel*(nLevel-1)/2;
  long long int fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::ivec samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::imat output(n,2);
  for(long long int i=0; i<n; ++i){
    output(i,0) = mapRow(samples(i),nLevel);
    output(i,1) = mapCol(samples(i),nLevel);
  }
  if(fullComb>0){
    arma::imat tmp(N*fullComb,2);
    long long int i;
    for(long long int j=0; j<(N*fullComb); ++j){
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
void changeId(Rcpp::IntegerVector newId,
              Rcpp::IntegerVector& oldId){
  oldId[0] = newId[0];
}

