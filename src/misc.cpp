// These functions may be called by R, but are not listed in the package namespace
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
//' @keywords internal
// [[Rcpp::export]]
arma::mat popVarCpp(const arma::mat& X) {
  if(X.n_rows==1){
    return(arma::mat(X.n_cols,X.n_cols,arma::fill::zeros));
  }else{
    return arma::cov(X,1);
  }
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

// The number of blocks used for nItem work items. Blocks are never empty,
// and there is always at least one, so that a function allocating one
// accumulator per block always has somewhere to put its results.
arma::uword countBlocks(arma::uword nItem){
  if(nItem<1){
    return 1;
  }
  if(nItem<nWorkBlocks){
    return nItem;
  }
  return nWorkBlocks;
}

// The first work item of a block. Called with block and block+1 to get the
// half open range a block covers. Any remainder is spread over the blocks
// instead of landing entirely on the last one.
arma::uword blockStart(arma::uword nItem, arma::uword nBlock,
                       arma::uword block){
  return (nItem*block)/nBlock;
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

// Samples random pairs without replacement from all possible combinations
// nLevel1 = number of levels for the first column
// nLevel2 = number of levels for the second column
// n = number of combinations to sample
// If n is larger than the total number of possible combinations (N), 
// then only n%N combinations are sampled and the rest are systematically assigned
// Returns an integer matrix with the sampled levels for each column
// Values in column 1 range from 1 to nLevel1
// Values in column 2 range from 1 to nLevel2
// Uses alphasimrRng::sampleInt(), which is seeded from R's RNG.
// Not safe for use inside OpenMP.
// [[Rcpp::export]]
arma::umat sampAllComb(arma::uword nLevel1, arma::uword nLevel2, 
                       arma::uword n){
  arma::uword N = nLevel1*nLevel2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  // Sample n combinations from the full set of size N
  arma::uvec samples = alphasimrRng::sampleInt(n,N);
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
// Uses alphasimrRng::sampleInt(), which is seeded from R's RNG.
// Not safe for use inside OpenMP.
// [[Rcpp::export]]
arma::umat sampHalfDialComb(arma::uword nLevel, arma::uword n){
  arma::uword N = nLevel*(nLevel-1)/2;
  arma::uword fullComb = 0;
  while(n>N){
    n -= N;
    ++fullComb;
  }
  // Sample n combinations from the full set of size N
  arma::uvec samples = alphasimrRng::sampleInt(n,N);
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

// n choose k recursive formula
double choose(double n, double k){ 
  if(k==0) return 1;
  return (n*choose(n-1,k-1))/k;
}

//' @title Check if OpenMP is available
//'
//' @description Checks if OpenMP is available
//'
//' @return logical
//'
//' @seealso \code{vignette("parallelization", package="AlphaSimR")}
//'  for setup details and \code{\link{getNumThreads}}.
//'
//' @examples
//' isOpenMPAvailable()
//' getNumThreads()
//'
//' @export
// [[Rcpp::export]]
bool isOpenMPAvailable(){
#ifdef _OPENMP
  return true;
#endif
  return false;
}

//' @title Number of available threads
//'
//' @description
//' Gets the number of available threads by calling the OpenMP function
//' \code{omp_get_max_threads()}
//'
//' @return integer
//'
//' @seealso \code{vignette("parallelization", package="AlphaSimR")}
//'  for setup details and \code{\link{isOpenMPAvailable}}.
//'
//' @examples
//' isOpenMPAvailable()
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
