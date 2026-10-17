#ifndef IBD_H
#define IBD_H

#include <RcppArmadillo.h>

// Recombination histories and IBD haplotypes are both stored per individual
// as a set of two-column integer matrices, one for each chromosome copy of
// each chromosome. The first column names a source haplotype and the second
// the site at which it starts. Storing these as nested R lists costs about
// thirty R objects per individual, which made the history about twenty times
// larger than its data and slowed every garbage collection. They are instead
// packed into a single integer vector per individual:
//
//   [0]               -nChr, negative so a packed record can be told apart
//                     from a founder's record, which is a vector of positive
//                     haplotype numbers
//   [1]               ploidy
//   [2, 2+nChr*ploidy)  number of rows for each chromosome copy, ordered by
//                     chromosome and then by copy
//   then              each matrix in the same order, stored column major,
//                     so its source column is followed by its site column
Rcpp::IntegerVector packRecord(
    const arma::field<arma::field<arma::Mat<int> > >& x);

#endif
