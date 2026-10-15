#ifndef PACKEDGV_H
#define PACKEDGV_H

// Genetic values worked out directly from packed genotypes. The functions
// are defined in getGv.cpp and shared with calcGenParam.cpp.

// The genotypes of a population, read in place from the R list holding them.
// Converting that list to an arma::field would copy every chromosome before
// any work was done, and decoding it to a dosage matrix, as getGeno does,
// would write out a byte for every individual and locus. Genetic values are
// worked out from the packed bits instead.
struct PackedGeno{
  std::vector<const unsigned char*> chrPtr;
  std::vector<arma::uword> nBins;
  arma::uword ploidy;
  arma::uword nInd;

  explicit PackedGeno(const Rcpp::S4& pop){
    Rcpp::List geno = pop.slot("geno");
    ploidy = Rcpp::as<arma::uword>(pop.slot("ploidy"));
    nInd = Rcpp::as<arma::uword>(pop.slot("nInd"));
    for(R_xlen_t i=0; i<geno.size(); ++i){
      // The pointers below are only valid while they point into the
      // population's own vectors, so a vector that would need coercing,
      // and would therefore be a temporary copy, is refused
      SEXP chrGeno = geno[i];
      if(TYPEOF(chrGeno)!=RAWSXP){
        Rcpp::stop("genotypes are not stored as raw vectors");
      }
      Rcpp::IntegerVector dim = Rf_getAttrib(chrGeno, R_DimSymbol);
      if(dim.size()!=3){
        Rcpp::stop("genotypes for a chromosome are not a three dimensional array");
      }
      chrPtr.push_back(RAW(chrGeno));
      nBins.push_back(arma::uword(dim[0]));
    }
  }
};

// Where each of a trait's loci sits in the packed genotypes: a pointer to its
// byte on the first haplotype of the first individual, the distance between
// haplotypes and between individuals, and the bit within the byte
struct LociPos{
  std::vector<const unsigned char*> base;
  arma::uvec hapStride;
  arma::uvec indStride;
  arma::uvec chr;
  arma::uvec byte;
  arma::uvec bit;
};

// Allele dosage of one individual at one locus
inline unsigned int dosageAt(const unsigned char* p, arma::uword hapStride,
                             arma::uword ploidy, arma::uword bit){
  unsigned int dose = 0;
  for(arma::uword h=0; h<ploidy; ++h){
    dose += (unsigned int) ((p[hapStride*h]>>bit) & 1u);
  }
  return dose;
}

// The effects of one set of genetic values on a set of loci, coded as the
// trait classes code them: additive effects on (x-P/2)*(2/P), dominance
// effects on x*(P-x)*(2/P)^2 and GxE effects on the additive coding, for a
// dosage x at ploidy P. A trait without epistasis or parent of origin effects
// is one of these, and so is any other function of dosage that has this form.
// hasA is false when there are no additive effects at all, so that their
// lookups can be skipped. Traits always have them.
struct GvEffects{
  arma::vec add;
  arma::vec dom;
  arma::vec gxe;
  bool hasA;
  bool hasD;
  bool hasGxe;
  double intercept;
  double gxeInt;
};

// Reads the dosages of a set of loci for a tile of individuals at a time.
// decode writes the dosage of locus i for the j-th individual of the tile
// to buf[i*tile+j], so that the dosages of one locus across the tile are
// contiguous. Loci sharing a byte are read together: spread[b] puts bit k of
// b into byte k of a 64 bit word, so summing it over the haplotypes gives
// the dosages of all eight loci of the byte at once.
class DosageTiles{
public:
  DosageTiles(const PackedGeno& g, const LociPos& pos);
  void decode(arma::uword indStart, arma::uword nTileInd, arma::uword tile,
              unsigned char* buf) const;
  arma::uword nLoci() const { return nLoci_; }
  // Number of distinct bytes the loci are in
  arma::uword nBytes() const { return keyBase_.size(); }
private:
  arma::uword ploidy_;
  arma::uword nLoci_;
  std::vector<uint64_t> spread_;
  std::vector<const unsigned char*> keyBase_;
  std::vector<arma::uword> keyHapStride_, keyIndStride_, keyStart_;
  std::vector<arma::uword> keyLocus_, keyShift_;
};

LociPos findLoci(const PackedGeno& g,
                 const arma::Col<int>& lociPerChr,
                 const arma::uvec& lociLoc);

// Genetic values for sets of effects that all have exactly the loci in pos.
// Each result holds the genetic values and, when there are GxE effects, the
// GxE values. Results do not depend on the thread count, and a set of
// effects gives the same values whichever others it is worked out with, as
// long as the same method is chosen for the group.
std::vector<arma::field<arma::vec> > gvGroupEffects(const std::vector<GvEffects>& effs,
                                                    const PackedGeno& g,
                                                    const LociPos& pos,
                                                    int nThreads);

#endif
