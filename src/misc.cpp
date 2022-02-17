// These functions may be called by R, but are not listed in the package namespace
#include "alphasimr.h"

std::bitset<8> toBits(unsigned char byte){
  return std::bitset<8>(byte);
}

unsigned char toByte(std::bitset<8> bits){
  return bits.to_ulong(); 
}

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
  arma::uword nChr = x.n_elem;
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
                                                      arma::uvec nBin,
                                                      arma::uword ploidy){
  arma::field<arma::Cube<unsigned char> > output(nBin.n_elem);
  arma::uword nTot = sum(nInd);
  arma::uword nPop = nInd.n_elem;
  // Allocate output
  for(arma::uword chr=0; chr<nBin.n_elem; ++chr){
    output(chr).set_size(nBin(chr),ploidy,nTot);
  }
  // Add individual genotypes
  arma::uword startInd=0, endInd=0;
  for(arma::uword i=0; i<nPop; ++i){
    if(nInd(i)>0){
      endInd += nInd(i)-1;
      Rcpp::S4 pop = popList[i];
      arma::field<arma::Cube<unsigned char> >geno = pop.slot("geno");
      for(arma::uword chr=0; chr<nBin.n_elem; ++chr){
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
arma::uword mapRow(const arma::uword& k, const arma::uword& n){
  return n-2-static_cast<arma::uword>(sqrt(-8*double(k) + 4*double(n)*(double(n)-1)-7)/2-0.5);
}

// Find column given mapping index
// row = previously determined row
// k = mapping index
// n = dimension of matrix
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n){
  return k+row+1 - n*(n-1)/2 + (n-row)*((n-row)-1)/2;
}


// Randomly samples integers without replacement
// n number of integers to return
// N number of integers to sample from
// Returns an integer vector of length n with values ranging from 0 to N-1
// Uses Jeffrey Scott Vitter's Method D
// [[Rcpp::export]]
arma::uvec sampleInt(arma::uword n, arma::uword N){
  arma::uvec output;
  output.set_size(n);
  if(n == 0){
    return output;
  }
  double q, v, x, y1, y2;
  arma::uword threshold = 13*n;
  arma::uword S, limit, top, bottom;
  arma::vec u(1,arma::fill::randu);
  v = exp(log(u(0))/double(n));
  q = double(N-n+1);
  while((n>1) & (threshold<N)){
    while(true){
      while(true){
        x = double(N)*(1-v);
        S = floor(x);
        if(double(S)<q){
          break;
        }
        u.randu();
        v = exp(log(u(0))/double(n));
      }
      u.randu();
      y1 = exp(log(u(0)*double(N)/q)/double(n-1));
      v = y1*(1-x/double(N))*(q/(q-double(S)));
      if(v <= 1){
        break;
      }
      y2 = 1;
      top = N-1;
      if((n-1) > S){
        bottom = N-n;
        limit = N-S;
      }else{
        bottom = N-S-1;
        limit = N-n+1;
      }
      for(arma::uword i=N-1; i>=limit; --i)
        y2 *= double(top)/double(bottom);
      u.randu();
      if((double(N)/(double(N)-x)) >= (y1*exp(log(y2)/double(n-1)))){
        v = exp(log(u(0))/double(n-1));
        break;
      }
      v = exp(log(u(0))/double(n));
    }
    output(n-1) = S+1;
    N = N-S-1;
    --n;
    q = double(N-n+1);
    threshold -= 13;
  }
  if(n > 1){
    top = N-n;
    while(n >= 2){
      u.randu();
      S = 0;
      q = double(top)/double(N);
      while(q > u(0)){
        ++S;
        --top;
        --N;
        q = (q*double(top))/double(N);
      }
      output(n-1) = S+1;
      --N;
      --n;
    }
    u.randu();
    output(0) = floor(u(0)*N);
  }else{
    output(0) = floor(v*N);
  }
  return cumsum(output);
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
arma::umat sampAllComb(arma::uword nLevel1, arma::uword nLevel2, 
                       arma::uword n){
  arma::uword N = nLevel1*nLevel2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::uvec samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::umat output(n,2);
  for(arma::uword  i=0; i<n; ++i){
    output(i,0) = samples(i)/nLevel2;
    output(i,1) = samples(i)%nLevel2;
  }
  if(fullComb>0){
    arma::umat tmp(N*fullComb,2);
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
arma::umat sampHalfDialComb(arma::uword nLevel, arma::uword n){
  arma::uword N = nLevel*(nLevel-1)/2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::uvec samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::umat output(n,2);
  for(arma::uword i=0; i<n; ++i){
    output(i,0) = mapRow(samples(i),nLevel);
    output(i,1) = mapCol(output(i,0),samples(i),nLevel);
  }
  if(fullComb>0){
    arma::umat tmp(N*fullComb,2);
    arma::uword i;
    for(arma::uword j=0; j<(N*fullComb); ++j){
      i = j%N;
      tmp(j,0) = mapRow(i,nLevel);
      tmp(j,1) = mapCol(tmp(j,0),i,nLevel);
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

// Knuth's algorithm for sampling from a Poisson distribution
arma::uword samplePoisson(double lambda){
  double p=1,L=exp(-lambda);
  arma::uword k=0;
  arma::vec u(1);
  do{
    k++;
    u.randu();
    p *= u(0);
  }while(p>L);
  return k-1;
}

// n choose k recursive formula
double choose(double n, double k){ 
  if(k==0) return 1;
  return (n*choose(n-1,k-1))/k;
}

// Gets the number of available threads
//' @title Number of available threads
//'
//' @description
//' Gets the number of available threads by calling the OpenMP function
//' \code{omp_get_max_threads()}
//'
//' @return integer
//'
//' @examples
//' getNumThreads()
//'
//' @export
// [[Rcpp::export]]
int getNumThreads(){
#ifdef _OPENMP
  return omp_get_max_threads();
#endif
  return 1;
}

