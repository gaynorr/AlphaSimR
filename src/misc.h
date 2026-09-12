#ifndef MISC_H
#define MISC_H

#include <bitset>

#include "rng.h"

// Splits a loop into a fixed number of work blocks.
//
// Several of AlphaSimR's parallel loops are run over blocks of work items
// rather than handing each thread a share of the items directly. The number
// of blocks is a fixed constant, so the split of the work, the number of
// accumulators a function allocates, and the order those accumulators are
// summed in, are the same on every machine. Results are then reproducible
// from set.seed() whatever the thread count happens to be.
const arma::uword nWorkBlocks = 64;

arma::uword countBlocks(arma::uword nItem);
arma::uword blockStart(arma::uword nItem, arma::uword nBlock,
                       arma::uword block);

arma::uword mapRow(const arma::uword& k, const arma::uword& n);
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n);
double choose(double n, double k);
// Defined here rather than in misc.cpp so that they can be inlined. They are
// called once per locus by the genotype accessors and once per byte by
// meiosis, and a call into another translation unit for each of those is a
// large share of the cost of both.
inline std::bitset<8> toBits(unsigned char byte){
  return std::bitset<8>(byte);
}

inline unsigned char toByte(std::bitset<8> bits){
  return (unsigned char) bits.to_ulong();
}

#endif
