#ifndef GETGENO_H
#define GETGENO_H

#include "packedGv.h"

// Allele dosage carried on haplotypes [pStart,pEnd) at a set of loci, read
// in place from packed genotypes. pStart=0 and pEnd=ploidy give the genotype
// dosage, and the first or second half of the haplotypes give the maternal
// or paternal dosage. lociLoc is one based.
arma::Mat<unsigned char> genoDosage(const PackedGeno& g,
                                    const arma::Col<int>& lociPerChr,
                                    arma::uvec lociLoc,
                                    arma::uword pStart, arma::uword pEnd,
                                    int nThreads);

// Haplotypes at a set of loci, as an nInd*ploidy by nLoci matrix
arma::Mat<unsigned char> genoHaplo(const PackedGeno& g,
                                   const arma::Col<int>& lociPerChr,
                                   arma::uvec lociLoc, int nThreads);

arma::mat genoToGenoA(const arma::Mat<unsigned char>& geno,
                      arma::uword ploidy, int nThreads);

arma::mat genoToGenoD(const arma::Mat<unsigned char>& geno,
                      arma::uword ploidy, int nThreads);

#endif
