#ifndef GENOVIEW_H
#define GENOVIEW_H

#include <vector>

// Packed genotypes seen as Armadillo cubes over memory that R owns.
//
// The meiosis functions read their parents' genotypes and write their
// progeny's a chromosome at a time through Cube::slice_colptr. Returning
// the progeny as an arma::field copied all of their genotypes on the way
// out, and taking the parents as one either copied the parents' genotypes
// or aliased them without marking them read only, as packedGv.h explains.
// The two classes below give the same cubes without a copy and with the
// parents' genotypes const.
//
// Each cube is built on memory it does not own (copy_aux_mem=false,
// strict=true), so it must never be resized, and the classes cannot be
// copied: copying an arma::Cube copies its elements into memory of its own,
// which would silently undo the point of the view. The vectors are reserved
// before any cube is added for the same reason, since growing them would
// copy the cubes too.

// Read only view of a list of raw arrays, such as pop@geno. The cubes alias
// R's memory, which other R objects may share, so they are only ever handed
// out as const. Nothing may be written through them.
class ConstGenoView{
public:
  explicit ConstGenoView(const Rcpp::List& geno){
    // PackedGeno checks that every chromosome is a raw array with three
    // dimensions and the same ploidy level and number of individuals
    PackedGeno g(geno);
    n_elem = g.chrPtr.size();
    cubes_.reserve(n_elem);
    for(arma::uword chr=0; chr<n_elem; ++chr){
      cubes_.emplace_back(const_cast<unsigned char*>(g.chrPtr[chr]),
                          g.nBins[chr], g.ploidy, g.nInd, false, true);
    }
  }
  ConstGenoView(const ConstGenoView&) = delete;
  ConstGenoView& operator=(const ConstGenoView&) = delete;

  const arma::Cube<unsigned char>& operator()(arma::uword chr) const{
    return cubes_[chr];
  }

  arma::uword n_elem;

private:
  std::vector<arma::Cube<unsigned char> > cubes_;
};

// Genotypes for new individuals, written in place into raw arrays that are
// allocated here and returned to R as they are. The arrays are allocated on
// the calling thread, because allocation calls R, and the cubes can then be
// written from any thread. Their contents start unset, as they did when the
// output was an arma::field sized with set_size, so every byte has to be
// written before list() is returned.
class NewGeno{
public:
  NewGeno(const std::vector<arma::uword>& nBins, arma::uword ploidy,
          arma::uword nInd){
    arma::uword nChr = nBins.size();
    list_ = Rcpp::List(nChr);
    cubes_.reserve(nChr);
    for(arma::uword chr=0; chr<nChr; ++chr){
      Rcpp::RawVector chrGeno = newChrGeno(nBins[chr], ploidy, nInd);
      list_[chr] = chrGeno;
      cubes_.emplace_back(RAW(chrGeno), nBins[chr], ploidy, nInd,
                          false, true);
    }
  }
  NewGeno(const NewGeno&) = delete;
  NewGeno& operator=(const NewGeno&) = delete;

  arma::Cube<unsigned char>& operator()(arma::uword chr){
    return cubes_[chr];
  }

  // The genotypes as an R list of raw arrays
  const Rcpp::List& list() const{
    return list_;
  }

private:
  Rcpp::List list_;
  std::vector<arma::Cube<unsigned char> > cubes_;
};

// Bytes per haplotype of each chromosome of a view, for sizing the output
inline std::vector<arma::uword> viewBins(const ConstGenoView& geno){
  std::vector<arma::uword> nBins(geno.n_elem);
  for(arma::uword chr=0; chr<geno.n_elem; ++chr){
    nBins[chr] = geno(chr).n_rows;
  }
  return nBins;
}

#endif
