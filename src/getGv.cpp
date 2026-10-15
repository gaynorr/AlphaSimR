#include "alphasimr.h"

// Calculates genetic values for genomic predictions using parental origin
arma::field<arma::vec> getGvA2(const Rcpp::S4& trait, 
                               const Rcpp::S4& pop, 
                               int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a1,a2,d;
  a1 = Rcpp::as<arma::vec>(trait.slot("addEff"));
  a2 = Rcpp::as<arma::vec>(trait.slot("addEffMale"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(a1.n_elem);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks);
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  output.set_size(1);
  output(0).set_size(nInd);
  // Half ploidy for xa
  arma::vec xa(ploidy/2+1);
  for(arma::uword i=0; i<xa.n_elem; ++i)
    xa(i) = (double(i)-dP/4.0)*(4.0/dP);
  // Full ploidy level for xd
  arma::vec xd(ploidy+1);
  for(arma::uword i=0; i<xd.n_elem; ++i)
    xd(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
  
  
  arma::Mat<unsigned char> maternalGeno = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> paternalGeno = getPaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(a1.n_elem, nBlocks, tid);
    arma::uword itemEnd = blockStart(a1.n_elem, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      arma::vec aEff1,aEff2,dEff;
      aEff1 = xa*a1(i);
      aEff2 = xa*a2(i);
      if(hasD){
        dEff = xd*d(i);
      }
      for(arma::uword j=0; j<nInd; ++j){
        gv(j,tid) += aEff1(maternalGeno(j,i)) + aEff2(paternalGeno(j,i));
        if(hasD){
          gv(j,tid) += dEff(maternalGeno(j,i)+paternalGeno(j,i));
        }
      }
    }
  }
  output(0) = sum(gv,1);
  return output;
}

// Individuals handled together in the inner loops below. A tile's genotypes
// stay in cache while every locus of a trait is swept past them. The value
// only affects speed: each individual sums its loci in the same order
// whatever the tile size, so results do not depend on it.
const arma::uword kGvIndTile = 32;

// Number of loci covered by one lookup table in gvByTable. Tables over four
// loci have sixteen entries and fit in two cache lines, so the tables for
// tens of thousands of loci stay in cache. Tables over a whole byte need one
// lookup per byte instead of two, but are sixteen times larger and were
// slower when several threads shared a cache.
const arma::uword kTabBits = 4;
const arma::uword kTabSize = arma::uword(1) << kTabBits;
const arma::uword kTabPerByte = 8/kTabBits;
const unsigned int kTabMask = (1u << kTabBits)-1u;

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

// Builds LociPos from the lociPerChr and lociLoc of a LociMap. Indexing is
// unchecked in the loops that use the result, so every locus is checked
// against the size of its chromosome here.
LociPos findLoci(const PackedGeno& g,
                 const arma::Col<int>& lociPerChr,
                 const arma::uvec& lociLoc){
  if(lociPerChr.n_elem!=g.chrPtr.size()){
    Rcpp::stop("number of chromosomes does not match the genotypes");
  }
  arma::uword nLoci = lociLoc.n_elem;
  if(arma::uword(arma::accu(lociPerChr))!=nLoci){
    Rcpp::stop("lociPerChr does not match the length of lociLoc");
  }
  LociPos pos;
  pos.base.resize(nLoci);
  pos.hapStride.set_size(nLoci);
  pos.indStride.set_size(nLoci);
  pos.chr.set_size(nLoci);
  pos.byte.set_size(nLoci);
  pos.bit.set_size(nLoci);
  arma::uword k = 0;
  for(arma::uword c=0; c<lociPerChr.n_elem; ++c){
    for(int j=0; j<lociPerChr(c); ++j, ++k){
      if(lociLoc(k)<1 || (lociLoc(k)-1)/8>=g.nBins[c]){
        Rcpp::stop("a locus lies outside its chromosome");
      }
      arma::uword loc = lociLoc(k)-1; // R to C++
      pos.chr(k) = c;
      pos.byte(k) = loc/8;
      pos.bit(k) = loc%8;
      pos.base[k] = g.chrPtr[c]+loc/8;
      pos.hapStride(k) = g.nBins[c];
      pos.indStride(k) = g.nBins[c]*g.ploidy;
    }
  }
  return pos;
}

// Picks out a trait's loci from the positions of all active QTL
LociPos subsetLoci(const LociPos& all, arma::Col<int> index){
  index -= 1; // R to C++
  LociPos pos;
  arma::uword n = index.n_elem;
  pos.base.resize(n);
  arma::uvec idx(n);
  for(arma::uword i=0; i<n; ++i){
    if(index(i)<0 || arma::uword(index(i))>=all.base.size()){
      Rcpp::stop("qtlIndex lies outside the active QTL");
    }
    idx(i) = arma::uword(index(i));
    pos.base[i] = all.base[idx(i)];
  }
  pos.hapStride = all.hapStride.elem(idx);
  pos.indStride = all.indStride.elem(idx);
  pos.chr = all.chr.elem(idx);
  pos.byte = all.byte.elem(idx);
  pos.bit = all.bit.elem(idx);
  return pos;
}

// Allele dosage of one individual at one locus
inline unsigned int dosageAt(const unsigned char* p, arma::uword hapStride,
                             arma::uword ploidy, arma::uword bit){
  unsigned int dose = 0;
  for(arma::uword h=0; h<ploidy; ++h){
    dose += (unsigned int) ((p[hapStride*h]>>bit) & 1u);
  }
  return dose;
}

// Sets up the output for a trait: genetic values, and GxE effects when the
// trait has them, each starting from its intercept
arma::field<arma::vec> startOutput(const Rcpp::S4& trait, arma::uword nInd,
                                   bool hasGxe, double gvShift,
                                   double gxeShift){
  arma::field<arma::vec> output(hasGxe ? 2 : 1);
  output(0).set_size(nInd);
  output(0).fill(double(trait.slot("intercept"))+gvShift);
  if(hasGxe){
    output(1).set_size(nInd);
    output(1).fill(double(trait.slot("gxeInt"))+gxeShift);
  }
  return output;
}

// Genetic values for additive traits, or additive and dominance traits in a
// diploid, using a lookup table for every group of kTabBits loci in the
// packed genotypes that holds one of the traits' loci.
//
// The additive coding (x-P/2)*(2/P) of a dosage x at ploidy P is linear in
// the alleles carried, so it splits into 2/P per allele on each haplotype
// less a constant. The constant is summed once and folded into the
// intercept, and a group of loci on one haplotype then contributes the table
// entry for its bit pattern. In a diploid the dominance coding is one exactly
// when the locus is heterozygous, which is when the bits of the two
// haplotypes differ, so the dominance table is read with their exclusive or.
// Above diploid the dominance coding is not a sum over haplotypes, and those
// traits go through gvByDosage instead.
//
// Every trait passed in must have exactly the loci in pos, as correlated
// traits added together do. They are worked out in one sweep, so the packed
// genotypes are read once for all of them, but each trait still sums its
// tables in the same order as it would alone, so grouping traits does not
// change their values. Tables are summed in a fixed order of chromosome and
// position, so results do not depend on the thread count, but they do differ
// from summing locus by locus in the last few bits.
std::vector<arma::field<arma::vec> > gvByTable(const std::vector<Rcpp::S4>& traits,
                                               const PackedGeno& g,
                                               const LociPos& pos,
                                               int nThreads){
  arma::uword nTrait = traits.size();
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  double alleleScale = 2.0/double(ploidy);
  arma::uword nLoci = pos.base.size();

  // One table per distinct group of loci, ordered by chromosome then
  // position. The key packs the two together so that sorting it gives that
  // order.
  arma::uvec maxByte(g.chrPtr.size()+1, arma::fill::zeros);
  for(arma::uword c=0; c<g.chrPtr.size(); ++c){
    maxByte(c+1) = maxByte(c)+g.nBins[c];
  }
  arma::uvec key(nLoci);
  for(arma::uword i=0; i<nLoci; ++i){
    key(i) = (maxByte(pos.chr(i))+pos.byte(i))*kTabPerByte+pos.bit(i)/kTabBits;
  }
  arma::uvec keys = arma::unique(key);
  arma::uword nKey = keys.n_elem;
  arma::uvec keyOf(nLoci);
  std::vector<const unsigned char*> keyBase(nKey);
  arma::uvec keyHapStride(nKey), keyIndStride(nKey);
  arma::Col<unsigned int> keyShift(nKey);
  for(arma::uword i=0; i<nLoci; ++i){
    arma::uword k = std::lower_bound(keys.begin(), keys.end(), key(i))-keys.begin();
    keyOf(i) = k;
    keyBase[k] = pos.base[i];
    keyHapStride(k) = pos.hapStride(i);
    keyIndStride(k) = pos.indStride(i);
    keyShift(k) = (unsigned int) ((pos.bit(i)/kTabBits)*kTabBits);
  }

  // Each table entry is the entry with its lowest set bit cleared plus that
  // bit's value, so a table costs one addition per entry
  arma::uvec lowBit(kTabSize, arma::fill::zeros);
  for(arma::uword p=1; p<kTabSize; ++p){
    arma::uword b = 0;
    while(((p>>b) & 1u)==0){
      ++b;
    }
    lowBit(p) = b;
  }

  // All of a key's tables, for every trait, sit next to each other in one
  // column, so a sweep reads them together. The offsets say where each
  // trait's tables start within a column, and are -1 for a table a trait
  // does not have.
  std::vector<bool> hasD(nTrait), hasGxe(nTrait);
  arma::ivec addOff(nTrait), domOff(nTrait), gxeOff(nTrait);
  arma::uword nTab = 0;
  for(arma::uword t=0; t<nTrait; ++t){
    hasD[t] = traits[t].hasSlot("domEff");
    hasGxe[t] = traits[t].hasSlot("gxeEff");
    addOff(t) = arma::sword(nTab++);
    domOff(t) = hasD[t] ? arma::sword(nTab++) : -1;
    gxeOff(t) = hasGxe[t] ? arma::sword(nTab++) : -1;
  }
  arma::mat tabs(kTabSize*nTab, nKey);
  std::vector<arma::field<arma::vec> > output(nTrait);
  arma::mat bitVal(kTabBits, nKey);
  auto buildTables = [&](arma::sword off){
    for(arma::uword k=0; k<nKey; ++k){
      double* tab = tabs.colptr(k)+kTabSize*arma::uword(off);
      const double* v = bitVal.colptr(k);
      tab[0] = 0.0;
      for(arma::uword p=1; p<kTabSize; ++p){
        tab[p] = tab[p & (p-1)]+v[lowBit(p)];
      }
    }
  };
  for(arma::uword t=0; t<nTrait; ++t){
    const Rcpp::S4& trait = traits[t];
    arma::vec a = Rcpp::as<arma::vec>(trait.slot("addEff"));
    if(a.n_elem!=nLoci){
      Rcpp::stop("number of loci does not match the number of effects");
    }
    // Per bit values first, so that two loci sharing a bit still both count
    double gvShift = 0.0, gxeShift = 0.0;
    bitVal.zeros();
    for(arma::uword i=0; i<nLoci; ++i){
      bitVal(pos.bit(i)%kTabBits,keyOf(i)) += a(i)*alleleScale;
      gvShift -= a(i);
    }
    buildTables(addOff(t));
    if(hasD[t]){
      arma::vec d = Rcpp::as<arma::vec>(trait.slot("domEff"));
      bitVal.zeros();
      for(arma::uword i=0; i<nLoci; ++i){
        bitVal(pos.bit(i)%kTabBits,keyOf(i)) += d(i);
      }
      buildTables(domOff(t));
    }
    if(hasGxe[t]){
      arma::vec ge = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
      bitVal.zeros();
      for(arma::uword i=0; i<nLoci; ++i){
        bitVal(pos.bit(i)%kTabBits,keyOf(i)) += ge(i)*alleleScale;
        gxeShift -= ge(i);
      }
      buildTables(gxeOff(t));
    }
    output[t] = startOutput(trait, nInd, hasGxe[t], gvShift, gxeShift);
  }
  std::vector<double*> gvPtr(nTrait), gxePtr(nTrait, NULL);
  for(arma::uword t=0; t<nTrait; ++t){
    gvPtr[t] = output[t](0).memptr();
    if(hasGxe[t]){
      gxePtr[t] = output[t](1).memptr();
    }
  }

  // Blocks of individuals own their part of the output, so the split is
  // fixed by the number of individuals and never by the thread count
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  if(nTrait==1){
    // A trait on its own is swept straight from the packed genotypes. This
    // loop is kept separate from the grouped one below because, with nothing
    // in it that varies by trait, it compiles to much faster code.
    double* gvOut = gvPtr[0];
    double* gxeOut = gxePtr[0];
    bool d0 = hasD[0], g0 = hasGxe[0];
    arma::uword aO = kTabSize*arma::uword(addOff(0));
    arma::uword dO = d0 ? kTabSize*arma::uword(domOff(0)) : 0;
    arma::uword gO = g0 ? kTabSize*arma::uword(gxeOff(0)) : 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword indStart = blockStart(nInd, nBlocks, block);
      arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
      for(arma::uword t0=indStart; t0<indEnd; t0+=kGvIndTile){
        arma::uword t1 = std::min(t0+kGvIndTile, indEnd);
        for(arma::uword k=0; k<nKey; ++k){
          const unsigned char* base = keyBase[k];
          arma::uword hs = keyHapStride(k);
          arma::uword is = keyIndStride(k);
          unsigned int sh = keyShift(k);
          const double* aTab = tabs.colptr(k)+aO;
          const double* dTab = tabs.colptr(k)+dO;
          const double* gTab = tabs.colptr(k)+gO;
          for(arma::uword j=t0; j<t1; ++j){
            const unsigned char* p = base+is*j;
            double s = 0.0;
            for(arma::uword h=0; h<ploidy; ++h){
              s += aTab[(p[hs*h]>>sh) & kTabMask];
            }
            if(d0){
              s += dTab[((p[0]^p[hs])>>sh) & kTabMask];
            }
            gvOut[j] += s;
            if(g0){
              double sg = 0.0;
              for(arma::uword h=0; h<ploidy; ++h){
                sg += gTab[(p[hs*h]>>sh) & kTabMask];
              }
              gxeOut[j] += sg;
            }
          }
        }
      }
    }
    return output;
  }

  // Traits sharing loci. For a tile of individuals and a chunk of keys, the
  // bit patterns are read out of the packed genotypes once into a buffer
  // small enough to stay in L1, and each trait then sweeps that buffer
  // contiguously. The larger tile means each trait's tables are read fewer
  // times. Neither size affects results: every trait still sums its keys in
  // order for each individual.
  const arma::uword kGroupTile = 256;
  const arma::uword kKeyChunk = 64;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    // Patterns by key, then haplotype, then individual, with the
    // exclusive or of a diploid's two haplotypes as an extra row
    arma::uword nRow = ploidy+(ploidy==2 ? 1 : 0);
    std::vector<unsigned char> buf(kKeyChunk*nRow*kGroupTile);
    for(arma::uword t0=indStart; t0<indEnd; t0+=kGroupTile){
      arma::uword t1 = std::min(t0+kGroupTile, indEnd);
      arma::uword nj = t1-t0;
      for(arma::uword k0=0; k0<nKey; k0+=kKeyChunk){
        arma::uword k1 = std::min(k0+kKeyChunk, nKey);
        for(arma::uword k=k0; k<k1; ++k){
          const unsigned char* base = keyBase[k];
          arma::uword hs = keyHapStride(k);
          arma::uword is = keyIndStride(k);
          unsigned int sh = keyShift(k);
          unsigned char* b = &buf[(k-k0)*nRow*kGroupTile];
          for(arma::uword j=0; j<nj; ++j){
            const unsigned char* p = base+is*(t0+j);
            for(arma::uword h=0; h<ploidy; ++h){
              b[h*kGroupTile+j] = (unsigned char) ((p[hs*h]>>sh) & kTabMask);
            }
            if(ploidy==2){
              b[2*kGroupTile+j] = (unsigned char) (((p[0]^p[hs])>>sh) & kTabMask);
            }
          }
        }
        for(arma::uword t=0; t<nTrait; ++t){
          double* gvOut = gvPtr[t]+t0;
          double* gxeOut = hasGxe[t] ? gxePtr[t]+t0 : NULL;
          bool dT = hasD[t], gT = hasGxe[t];
          arma::uword aO = kTabSize*arma::uword(addOff(t));
          arma::uword dO = dT ? kTabSize*arma::uword(domOff(t)) : 0;
          arma::uword gO = gT ? kTabSize*arma::uword(gxeOff(t)) : 0;
          for(arma::uword k=k0; k<k1; ++k){
            const double* aTab = tabs.colptr(k)+aO;
            const double* dTab = tabs.colptr(k)+dO;
            const double* gTab = tabs.colptr(k)+gO;
            const unsigned char* b = &buf[(k-k0)*nRow*kGroupTile];
            for(arma::uword j=0; j<nj; ++j){
              double s = 0.0;
              for(arma::uword h=0; h<ploidy; ++h){
                s += aTab[b[h*kGroupTile+j]];
              }
              if(dT){
                s += dTab[b[2*kGroupTile+j]];
              }
              gvOut[j] += s;
              if(gT){
                double sg = 0.0;
                for(arma::uword h=0; h<ploidy; ++h){
                  sg += gTab[b[h*kGroupTile+j]];
                }
                gxeOut[j] += sg;
              }
            }
          }
        }
      }
    }
  }
  return output;
}

// Genetic values for additive or additive and dominance traits at any
// ploidy, looking up each locus's value from its dosage.
//
// Every trait passed in must have exactly the loci in pos. For a tile of
// individuals and a chunk of loci, dosages are read out of the packed
// genotypes once into a buffer small enough to stay in L1, and each trait
// then sweeps that buffer contiguously with a table over dosages for each
// locus. Each individual sums its loci in trait order, which is the order
// the dosage matrix version this replaces used, so results match it exactly.
std::vector<arma::field<arma::vec> > gvByDosage(const std::vector<Rcpp::S4>& traits,
                                                const PackedGeno& g,
                                                const LociPos& pos,
                                                int nThreads){
  arma::uword nTrait = traits.size();
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  arma::uword nGeno = ploidy+1;
  double dP = double(ploidy);
  arma::uword nLoci = pos.base.size();
  arma::vec x(nGeno); // Genotype dosage
  for(arma::uword i=0; i<nGeno; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);

  // What a locus contributes depends only on its dosage, so the value of
  // each dosage is worked out once per locus and looked up afterward
  std::vector<arma::mat> effTab(nTrait), gxeTab(nTrait);
  std::vector<bool> hasGxe(nTrait);
  std::vector<arma::field<arma::vec> > output(nTrait);
  for(arma::uword t=0; t<nTrait; ++t){
    const Rcpp::S4& trait = traits[t];
    bool hasD = trait.hasSlot("domEff");
    hasGxe[t] = trait.hasSlot("gxeEff");
    arma::vec a = Rcpp::as<arma::vec>(trait.slot("addEff"));
    if(a.n_elem!=nLoci){
      Rcpp::stop("number of loci does not match the number of effects");
    }
    arma::vec d;
    if(hasD){
      d = Rcpp::as<arma::vec>(trait.slot("domEff"));
    }
    effTab[t].set_size(nGeno,nLoci);
    for(arma::uword i=0; i<nLoci; ++i){
      if(hasD){
        effTab[t].col(i) = xa*a(i)+xd*d(i);
      }else{
        effTab[t].col(i) = xa*a(i);
      }
    }
    if(hasGxe[t]){
      arma::vec ge = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
      gxeTab[t].set_size(nGeno,nLoci);
      for(arma::uword i=0; i<nLoci; ++i){
        gxeTab[t].col(i) = xa*ge(i);
      }
    }
    output[t] = startOutput(trait, nInd, hasGxe[t], 0.0, 0.0);
  }

  // Neither size affects results, only speed
  const arma::uword kTile = 256;
  const arma::uword kChunk = 64;
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    std::vector<unsigned char> buf(kChunk*kTile);
    for(arma::uword t0=indStart; t0<indEnd; t0+=kTile){
      arma::uword nj = std::min(t0+kTile, indEnd)-t0;
      for(arma::uword i0=0; i0<nLoci; i0+=kChunk){
        arma::uword i1 = std::min(i0+kChunk, nLoci);
        for(arma::uword i=i0; i<i1; ++i){
          const unsigned char* base = pos.base[i];
          arma::uword hs = pos.hapStride(i);
          arma::uword is = pos.indStride(i);
          arma::uword bit = pos.bit(i);
          unsigned char* b = &buf[(i-i0)*kTile];
          for(arma::uword j=0; j<nj; ++j){
            b[j] = (unsigned char) dosageAt(base+is*(t0+j), hs, ploidy, bit);
          }
        }
        for(arma::uword t=0; t<nTrait; ++t){
          double* gvOut = output[t](0).memptr()+t0;
          for(arma::uword i=i0; i<i1; ++i){
            const double* eff = effTab[t].colptr(i);
            const unsigned char* b = &buf[(i-i0)*kTile];
            for(arma::uword j=0; j<nj; ++j){
              gvOut[j] += eff[b[j]];
            }
          }
          if(hasGxe[t]){
            double* gxeOut = output[t](1).memptr()+t0;
            for(arma::uword i=i0; i<i1; ++i){
              const double* gEff = gxeTab[t].colptr(i);
              const unsigned char* b = &buf[(i-i0)*kTile];
              for(arma::uword j=0; j<nj; ++j){
                gxeOut[j] += gEff[b[j]];
              }
            }
          }
        }
      }
    }
  }
  return output;
}

// Genetic values for a trait with additive-by-additive epistasis. Every
// locus belongs to exactly one pair, and a pair's additive, dominance and
// epistatic terms depend only on its two dosages, so each pair has a table
// over the two dosages. Each individual sums its pairs in order, as the
// dosage matrix version this replaces did, so results match it exactly.
arma::field<arma::vec> gvEpistasis(const Rcpp::S4& trait, const PackedGeno& g,
                                   const LociPos& pos, int nThreads){
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  arma::uword nGeno = ploidy+1;
  double dP = double(ploidy);
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
  arma::vec a,d,ge;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  if(hasGxe){
    ge = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
  }
  if(pos.base.size()!=a.n_elem){
    Rcpp::stop("number of loci does not match the number of effects");
  }
  arma::uword nPair = E.n_rows;
  arma::vec x(nGeno); // Genotype dosage
  for(arma::uword i=0; i<nGeno; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);

  arma::mat effTab(nGeno*nGeno,nPair), gxeTab;
  if(hasGxe){
    gxeTab.set_size(nGeno*nGeno,nPair);
  }
  arma::uvec loc1(nPair), loc2(nPair);
  for(arma::uword i=0; i<nPair; ++i){
    if(E(i,0)<0 || E(i,1)<0 || E(i,0)>=a.n_elem || E(i,1)>=a.n_elem){
      Rcpp::stop("epiEff refers to a locus the trait does not have");
    }
    arma::uword l1 = arma::uword(E(i,0));
    arma::uword l2 = arma::uword(E(i,1));
    loc1(i) = l1;
    loc2(i) = l2;
    double a1 = a(l1), a2 = a(l2), e12 = E(i,2);
    double d1 = 0.0, d2 = 0.0, g1 = 0.0, g2 = 0.0;
    if(hasD){
      d1 = d(l1);
      d2 = d(l2);
    }
    if(hasGxe){
      g1 = ge(l1);
      g2 = ge(l2);
    }
    double* eff = effTab.colptr(i);
    double* gEff = hasGxe ? gxeTab.colptr(i) : NULL;
    for(arma::uword u=0; u<nGeno; ++u){
      for(arma::uword v=0; v<nGeno; ++v){
        if(hasD){
          eff[u+v*nGeno] = a1*xa(u)+d1*xd(u)+a2*xa(v)+d2*xd(v)+e12*xa(u)*xa(v);
        }else{
          eff[u+v*nGeno] = a1*xa(u)+a2*xa(v)+e12*xa(u)*xa(v);
        }
        if(hasGxe){
          gEff[u+v*nGeno] = g1*xa(u)+g2*xa(v);
        }
      }
    }
  }

  arma::field<arma::vec> output = startOutput(trait, nInd, hasGxe, 0.0, 0.0);
  double* gvPtr = output(0).memptr();
  double* gxePtr = hasGxe ? output(1).memptr() : NULL;

  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword t0=indStart; t0<indEnd; t0+=kGvIndTile){
      arma::uword t1 = std::min(t0+kGvIndTile, indEnd);
      for(arma::uword i=0; i<nPair; ++i){
        arma::uword l1 = loc1(i), l2 = loc2(i);
        const unsigned char* base1 = pos.base[l1];
        const unsigned char* base2 = pos.base[l2];
        arma::uword hs1 = pos.hapStride(l1), hs2 = pos.hapStride(l2);
        arma::uword is1 = pos.indStride(l1), is2 = pos.indStride(l2);
        arma::uword bit1 = pos.bit(l1), bit2 = pos.bit(l2);
        const double* eff = effTab.colptr(i);
        const double* gEff = hasGxe ? gxeTab.colptr(i) : NULL;
        for(arma::uword j=t0; j<t1; ++j){
          unsigned int u = dosageAt(base1+is1*j, hs1, ploidy, bit1);
          unsigned int v = dosageAt(base2+is2*j, hs2, ploidy, bit2);
          gvPtr[j] += eff[u+v*nGeno];
          if(hasGxe){
            gxePtr[j] += gEff[u+v*nGeno];
          }
        }
      }
    }
  }
  return output;
}

// Works out genetic values for traits that share exactly the same loci and
// have no epistasis or parent of origin effects, choosing between the two
// methods above.
//
// gvByTable does one lookup per haplotype per table for each trait, plus one
// for dominance and one per haplotype for GxE, and a table covers kTabBits
// loci. gvByDosage reads each locus's dosage once, which costs a little more
// than a table lookup per haplotype, and then does one cheaper lookup per
// locus for each trait, plus one for GxE. Measured on a diploid, a table
// lookup costs about c, reading a dosage 1.25c per haplotype and a dosage
// lookup 0.5c, and the choice below compares the two totals. Tables win for
// a single trait whatever the density, because the dosages are read for that
// one trait alone, while for many traits sharing sparse loci the cheaper
// dosage lookups win. Polyploid dominance can only use dosages. The choice
// depends only on the traits and their loci, never on the thread count.
std::vector<arma::field<arma::vec> > gvGroup(const std::vector<Rcpp::S4>& traits,
                                             const PackedGeno& g,
                                             const LociPos& pos,
                                             int nThreads){
  double P = double(g.ploidy);
  bool needDosage = false;
  double tableWork = 0.0, dosageWork = 1.25*P;
  for(arma::uword t=0; t<traits.size(); ++t){
    bool hasD = traits[t].hasSlot("domEff");
    bool hasGxe = traits[t].hasSlot("gxeEff");
    if(hasD && g.ploidy!=2){
      needDosage = true;
    }
    tableWork += P+(hasD ? 1.0 : 0.0)+(hasGxe ? P : 0.0);
    dosageWork += 0.5*(1.0+(hasGxe ? 1.0 : 0.0));
  }
  arma::uword nLoci = pos.base.size();
  if(!needDosage && nLoci>0){
    // Number of distinct tables the loci would need
    std::vector<std::pair<arma::uword,arma::uword> > keys(nLoci);
    for(arma::uword i=0; i<nLoci; ++i){
      keys[i] = std::make_pair(pos.chr(i), pos.byte(i)*kTabPerByte+pos.bit(i)/kTabBits);
    }
    std::sort(keys.begin(), keys.end());
    arma::uword nKey = std::unique(keys.begin(), keys.end())-keys.begin();
    bool useTable = double(nKey)*tableWork < double(nLoci)*dosageWork;
    if(useTable){
      return gvByTable(traits, g, pos, nThreads);
    }
  }
  return gvByDosage(traits, g, pos, nThreads);
}

// Chooses how to work out a trait's genetic values once its loci are found.
// Traits with parent of origin effects are handled by the callers.
arma::field<arma::vec> gvFromPacked(const Rcpp::S4& trait, const PackedGeno& g,
                                    const LociPos& pos, int nThreads){
  if(trait.hasSlot("epiEff")){
    return gvEpistasis(trait, g, pos, nThreads);
  }
  return gvGroup(std::vector<Rcpp::S4>(1, trait), g, pos, nThreads)[0];
}
// Calculates genetic values for a trait
// Returns output in a list with length 1 or 2
//   The first item contains genetic values
//   The second item contains GxE effects (optional)
// [[Rcpp::export]]
arma::field<arma::vec> getGv(const Rcpp::S4& trait,
                             const Rcpp::S4& pop,
                             int nThreads){
  if(trait.hasSlot("addEffMale")){
    // Genomic prediction
    return getGvA2(trait, pop, nThreads);
  }
  PackedGeno g(pop);
  LociPos pos = findLoci(g, Rcpp::as<arma::Col<int> >(trait.slot("lociPerChr")),
                         Rcpp::as<arma::uvec>(trait.slot("lociLoc")));
  return gvFromPacked(trait, g, pos, nThreads);
}

// Calculates genetic values for genomic predictions using parental origin
arma::field<arma::vec> getGvIndexA2(const Rcpp::S4& trait, 
                                    const Rcpp::S4& pop, 
                                    int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a1,a2,d;
  a1 = Rcpp::as<arma::vec>(trait.slot("addEff"));
  a2 = Rcpp::as<arma::vec>(trait.slot("addEffMale"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(a1.n_elem);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks);
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  output.set_size(1);
  output(0).set_size(nInd);
  // Half ploidy for xa
  arma::vec xa(ploidy/2+1);
  for(arma::uword i=0; i<xa.n_elem; ++i)
    xa(i) = (double(i)-dP/4.0)*(4.0/dP);
  // Full ploidy level for xd
  arma::vec xd(ploidy+1);
  for(arma::uword i=0; i<xd.n_elem; ++i)
    xd(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
  
  
  arma::Mat<unsigned char> maternalGeno = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> paternalGeno = getPaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(a1.n_elem, nBlocks, tid);
    arma::uword itemEnd = blockStart(a1.n_elem, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      arma::vec aEff1,aEff2,dEff;
      aEff1 = xa*a1(i);
      aEff2 = xa*a2(i);
      if(hasD){
        dEff = xd*d(i);
      }
      for(arma::uword j=0; j<nInd; ++j){
        gv(j,tid) += aEff1(maternalGeno(j,i)) + aEff2(paternalGeno(j,i));
        if(hasD){
          gv(j,tid) += dEff(maternalGeno(j,i)+paternalGeno(j,i));
        }
      }
    }
  }
  output(0) = sum(gv,1);
  return output;
}

// [[Rcpp::export]]
Rcpp::List getGvIndex(const Rcpp::S4& pop,
                      const Rcpp::List& traitList,
                      const Rcpp::S4& activeQtl,
                      const arma::field<arma::Col<int> >& qtlIndex,
                      arma::uword nTraits,
                      int nThreads){

  // Get summary variables
  arma::uword nInd = pop.slot("nInd");

  // The active QTL are located once and each trait picks out its own loci.
  // No genotypes are decoded here: each trait reads the packed genotypes
  // for just the loci it uses.
  PackedGeno g(pop);
  LociPos allPos = findLoci(g, Rcpp::as<arma::Col<int> >(activeQtl.slot("lociPerChr")),
                            Rcpp::as<arma::uvec>(activeQtl.slot("lociLoc")));

  // Format output
  arma::mat gv(nInd, nTraits);
  arma::field<arma::vec> gxe(nTraits);


  // Traits without epistasis that have exactly the same loci, as correlated
  // traits added together do, are grouped so that one sweep of the
  // genotypes serves all of them. Grouping does not change any trait's
  // values, only how often the genotypes are read.
  std::vector<arma::field<arma::vec> > outputs(nTraits);
  std::vector<std::vector<arma::uword> > groups;
  for(arma::uword i=0; i<nTraits; ++i){
    const Rcpp::S4& trait = traitList[i];
    if(trait.hasSlot("addEffMale")){
      // Trait with separate male and female effects (from genomic prediction)
      // Unlikely to be use
      outputs[i] = getGvIndexA2(trait, pop, nThreads);
    }else if(trait.hasSlot("epiEff")){
      outputs[i] = gvFromPacked(trait, g, subsetLoci(allPos, qtlIndex(i)),
                                nThreads);
    }else{
      bool placed = false;
      for(arma::uword k=0; k<groups.size() && !placed; ++k){
        const arma::Col<int>& other = qtlIndex(groups[k][0]);
        if(other.n_elem==qtlIndex(i).n_elem &&
           arma::all(other==qtlIndex(i))){
          groups[k].push_back(i);
          placed = true;
        }
      }
      if(!placed){
        groups.push_back(std::vector<arma::uword>(1, i));
      }
    }
  }
  for(arma::uword k=0; k<groups.size(); ++k){
    std::vector<Rcpp::S4> members;
    for(arma::uword t=0; t<groups[k].size(); ++t){
      members.push_back(traitList[groups[k][t]]);
    }
    std::vector<arma::field<arma::vec> > groupOut =
      gvGroup(members, g, subsetLoci(allPos, qtlIndex(groups[k][0])),
              nThreads);
    for(arma::uword t=0; t<groups[k].size(); ++t){
      outputs[groups[k][t]] = groupOut[t];
    }
  }

  for(arma::uword i=0; i<nTraits; ++i){
    gv.col(i) = outputs[i](0);
    if(outputs[i].n_elem == 2){
      // Has GxE
      gxe(i) = outputs[i](1);
    }
  }

  return Rcpp::List::create(Rcpp::Named("gv") = gv,
                            Rcpp::Named("gxe") = gxe);
}
