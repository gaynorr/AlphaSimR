// Fast mode of the multikernel RR-BLUP models: RRBLUP_D, RRBLUP_GCA and
// RRBLUP_SCA with fastMode TRUE.
//
// These models have several random terms, each a matrix of marker codes with
// its own variance, or for several traits its own covariance matrix among the
// traits:
//
//   RRBLUP_D    additive dosage and dominance
//   RRBLUP_GCA  maternal and paternal haplotype dosage
//   RRBLUP_SCA  maternal and paternal haplotype dosage and dominance
//
// The fast mode of RRBLUP solves a single term, and for several traits it
// separates them by a canonical transformation of Ve and Vu. Neither carries
// over: every term has its own variance ratio, and two or three covariance
// matrices cannot be diagonalized together with Ve. The model is therefore
// solved here as one system over the effects of every term and trait, by
// preconditioned conjugate gradient, iterating on the genotypes rather than
// on a stored coefficient matrix as the fast mode of RRBLUP does (Stranden
// and Lidauer, 1999, doi:10.3168/jds.S0022-0302(99)75535-9).
//
// The system is written in whitened effects. The responses are multiplied by
// Le^-T, for Ve = Le*Le', which leaves errors with identity covariance, and
// each term's covariance on that scale, Le^-1*Vu_k*Le^-T, is written as
// F_k*F_k' over its non-negative eigenvalues. With U_k = A_k*F_k'*Le' the
// effects A_k are independent with unit variance, and the coefficient matrix
// is the identity plus F'*W'*S*W*F, which is positive definite whatever Vu_k
// is. A term or direction with no variance has no columns in F_k and drops
// out, with no special case, and a single trait is one column.
//
// The genotypes are read from the copy packed to one bit per allele that the
// fast mode of RRBLUP iterates on, and every term is read from it through
// lookup tables over four markers at a time. A dominance code, x*(P-x)
// scaled, counts the unordered pairs of haplotypes carrying different
// alleles, so it is a sum over pairs of haplotypes of the exclusive or of
// their bits, and is looked up the same way as an allele count.
//
// The preconditioner inverts each marker's block of the coefficient matrix,
// over every term and whitened trait, and deflates the largest eigenvalues
// that the relationships between individuals give the equations, as the
// fast mode of RRBLUP does; see mkDeflation.

#include <algorithm>
#include <utility>

#include "mme.h"

// Defined in mmeMultiKernel.cpp. Declared here without its default, which
// Rcpp reads from the definition.
Rcpp::List solveRRBLUPMK(const arma::mat& y, const arma::mat& X,
                         const arma::field<arma::mat>& Mlist,
                         int maxIter);
Rcpp::List remlMultiKernelGram(const arma::mat& Y, const arma::mat& X,
                               arma::field<arma::mat>& Glist,
                               const arma::vec& rank, int maxIter,
                               const char* name);
// Defined in MME.cpp, with the deflation of the fast mode of RRBLUP
arma::uvec deflationSubset(arma::uword nInd, arma::uword nSub);
void applyDeflation(arma::vec& v, const arma::mat& W, const arma::vec& shrink,
                    int nThreads);

namespace {

// Iterations allowed to the REML fit of the variance components on the
// subset of records, the default of maxIter in the direct form of these
// models
const int fastMKRemlIter = 100;

enum class TermKind { additive, maternal, paternal, dominance };

// How a random term is read from the packed genotypes. An additive term sums
// the alleles on haplotypes [hapStart,hapEnd), and its code is
// scale*count-constant. A dominance term counts the pairs of haplotypes with
// different alleles, and its code is scale*count.
struct FastTerm{
  TermKind kind;
  arma::uword hapStart;
  arma::uword hapEnd;
  double scale;
  double constant;
};

// The codes match genoToGenoA and genoToGenoD, which the direct form of the
// models uses: (x-P/2)*(2/P) for a dosage x at ploidy P, the same with P/2
// for the dosage of either parent's haplotypes, and x*(P-x)*(2/P)^2 for
// dominance
FastTerm makeTerm(TermKind kind, arma::uword ploidy){
  FastTerm t;
  t.kind = kind;
  double dP = double(ploidy);
  arma::uword half = ploidy/2;
  switch(kind){
  case TermKind::additive:
    t.hapStart = 0;
    t.hapEnd = ploidy;
    t.scale = 2.0/dP;
    t.constant = 1.0;
    break;
  case TermKind::maternal:
    t.hapStart = 0;
    t.hapEnd = half;
    t.scale = 2.0/double(half);
    t.constant = 1.0;
    break;
  case TermKind::paternal:
    t.hapStart = half;
    t.hapEnd = ploidy;
    t.scale = 2.0/double(half);
    t.constant = 1.0;
    break;
  case TermKind::dominance:
    t.hapStart = 0;
    t.hapEnd = ploidy;
    t.scale = (2.0/dP)*(2.0/dP);
    t.constant = 0.0;
    break;
  }
  return t;
}

// The alleles of individual i at marker j, one per haplotype
inline void readBits(const arma::Mat<unsigned char>& planes,
                     arma::uword ploidy, arma::uword i, arma::uword j,
                     unsigned int* b){
  const unsigned char* col = planes.colptr(j/8);
  unsigned int s = (unsigned int) (j%8);
  for(arma::uword h=0; h<ploidy; ++h){
    b[h] = ((unsigned int) col[i*ploidy+h]>>s) & 1u;
  }
}

// A term's code for one individual at one marker, from its alleles
inline double termValue(const FastTerm& t, const unsigned int* b,
                        arma::uword ploidy){
  if(t.kind==TermKind::dominance){
    unsigned int x = 0;
    for(arma::uword h=0; h<ploidy; ++h){
      x += b[h];
    }
    // The number of unordered pairs of haplotypes with different alleles
    return t.scale*double(x*((unsigned int) ploidy-x));
  }
  unsigned int x = 0;
  for(arma::uword h=t.hapStart; h<t.hapEnd; ++h){
    x += b[h];
  }
  return t.scale*double(x)-t.constant;
}

// Builds the lookup tables of one term for the products below, for every
// column of the packed genotypes and every trait, as multMcBlock does. Entry
// q of the low half byte's table starts at q*nT and of the high one's at
// (16+q)*nT, from tab.
void buildTables(double* tabBase, arma::uword tabStride,
                 const arma::mat& Ut, double scale, arma::uword nCol){
  arma::uword nT = Ut.n_rows;
  arma::uword nLoci = Ut.n_cols;
  std::vector<double> v(8*nT);
  for(arma::uword c=0; c<nCol; ++c){
    for(arma::uword b=0; b<8; ++b){
      arma::uword j = 8*c+b;
      for(arma::uword k=0; k<nT; ++k){
        v[b*nT+k] = (j<nLoci) ? scale*Ut(k,j) : 0.0;
      }
    }
    double* tab = tabBase+c*tabStride;
    for(arma::uword half=0; half<2; ++half){
      double* tb = tab+16*half*nT;
      const double* vb = v.data()+4*half*nT;
      for(arma::uword k=0; k<nT; ++k){
        tb[k] = 0.0;
      }
      for(unsigned int q=1; q<16; ++q){
        unsigned int low = q & (0u-q);
        unsigned int bit = (low==1u) ? 0 : ((low==2u) ? 1 : ((low==4u) ? 2 : 3));
        const double* prev = tb+(q & (q-1u))*nT;
        const double* add = vb+bit*nT;
        double* out = tb+q*nT;
        for(arma::uword k=0; k<nT; ++k){
          out[k] = prev[k]+add[k];
        }
      }
    }
  }
}

// Adds one byte's contribution from a table to an individual's results
inline void addLookup(double* ti, const double* tab, unsigned int b,
                      arma::uword nT){
  const double* lo = tab+(b & 15u)*nT;
  const double* hi = tab+(16+(b>>4))*nT;
  for(arma::uword k=0; k<nT; ++k){
    ti[k] += lo[k]+hi[k];
  }
}

// The record loops of multTerms and multTermsT with the ploidy P, the model
// and the number of traits NT fixed at compile time. Every lookup's term and
// byte are then constants: a byte is a haplotype's, row[a], or for a
// dominance term the exclusive or of a pair, row[a]^row[b], read straight
// from the record's packed genotypes. An earlier form built a record's bytes
// in an array, by loops over the haplotypes and their pairs, and looked them
// up through a list made at run time. GCC at -O2 does not unroll those loops
// for a tetraploid, and the array they wrote to the stack, read back to
// address the bins, made the transpose two to three times slower than the
// general loop. Nothing here depends on the compiler unrolling a loop.
//
// The models are 0 for RRBLUP_D (additive and dominance), 1 for RRBLUP_GCA
// (maternal and paternal) and 2 for RRBLUP_SCA (maternal, paternal and
// dominance). Their lookups are made in the order the general loops make
// them: for each term in turn, an additive term looks up each of its
// haplotypes and a dominance term each pair of haplotypes, a before b. The
// additions are therefore the same, in the same order, and give the same
// result.

// One lookup: the term whose table and bins it uses, and its byte
struct FixedLook{
  std::size_t term;
  std::size_t a;
  std::size_t b;
  bool pair;
};

// The lookups of model M at ploidy P, in order
template<std::size_t P, int M>
struct FixedLooks{
  static constexpr std::size_t nPair = P*(P-1)/2;
  static constexpr std::size_t N = (M==1) ? P : P+nPair;
  struct Arr{
    FixedLook v[N];
  };
  static constexpr Arr make(){
    Arr r{};
    std::size_t n = 0;
    if(M==0){
      for(std::size_t h=0; h<P; ++h){
        r.v[n++] = {0, h, h, false};
      }
    }else{
      for(std::size_t h=0; h<P/2; ++h){
        r.v[n++] = {0, h, h, false};
      }
      for(std::size_t h=P/2; h<P; ++h){
        r.v[n++] = {1, h, h, false};
      }
    }
    if(M!=1){
      std::size_t t = (M==0) ? 1 : 2;
      for(std::size_t h=0; h+1<P; ++h){
        for(std::size_t h2=h+1; h2<P; ++h2){
          r.v[n++] = {t, h, h2, true};
        }
      }
    }
    return r;
  }
  static constexpr Arr looks = make();
};

// The byte lookup L of model M reads from a record's haplotypes
template<std::size_t P, int M, std::size_t L>
inline unsigned int fixedByte(const unsigned char* row){
  constexpr FixedLook l = FixedLooks<P,M>::looks.v[L];
  if constexpr(l.pair){
    return (unsigned int) (row[l.a]^row[l.b]);
  }else{
    return (unsigned int) row[l.a];
  }
}

// One record of multTermsFixed at one column. As in multMcRecord in
// MME.cpp, the loops over traits are pack expansions, so that every index
// into the accumulators is a constant and they can be kept in registers.
template<std::size_t P, int M, std::size_t... L, std::size_t... K>
inline void multTermsRecord(double* ti, const unsigned char* row,
                            const double* tabCol, std::index_sequence<L...>,
                            std::index_sequence<K...>){
  constexpr std::size_t NT = sizeof...(K);
  constexpr arma::uword termLen = 32*NT;
  double s[NT] = {ti[K]...};
  auto one = [&](auto lc){
    constexpr std::size_t l = decltype(lc)::value;
    unsigned int b = fixedByte<P,M,l>(row);
    const double* tab = tabCol+FixedLooks<P,M>::looks.v[l].term*termLen;
    const double* lo = tab+(b & 15u)*NT;
    const double* hi = tab+(16+(b>>4))*NT;
    ((s[K] += lo[K]+hi[K]), ...);
  };
  (one(std::integral_constant<std::size_t,L>()), ...);
  ((ti[K] = s[K]), ...);
}

// Records [indStart,indEnd) of multTerms, Tt holding NT values per record
template<std::size_t P, int M, std::size_t NT>
void multTermsFixed(double* tPtr, const arma::Mat<unsigned char>& planes,
                    const double* tabsPtr, arma::uword tabStride,
                    arma::uword indStart, arma::uword indEnd){
  const arma::uword nCol = planes.n_cols;
  for(arma::uword c=0; c<nCol; ++c){
    const unsigned char* col = planes.colptr(c);
    const double* tabCol = tabsPtr+c*tabStride;
    for(arma::uword i=indStart; i<indEnd; ++i){
      multTermsRecord<P,M>(tPtr+i*NT, col+i*P, tabCol,
                           std::make_index_sequence<FixedLooks<P,M>::N>{},
                           std::make_index_sequence<NT>{});
    }
  }
}

// One record of multTermsT at one column, into each term's two sets of bins
// of setLen = 32*NT values, the set of the record's parity
template<std::size_t P, int M, std::size_t NT, std::size_t... L>
inline void multTermsTRecord(double* binPtr, const unsigned char* row,
                             const double* ti, arma::uword parity,
                             std::index_sequence<L...>){
  constexpr arma::uword setLen = 32*NT;
  auto one = [&](auto lc){
    constexpr std::size_t l = decltype(lc)::value;
    unsigned int b = fixedByte<P,M,l>(row);
    double* bt = binPtr+(2*FixedLooks<P,M>::looks.v[l].term+parity)*setLen;
    double* lo = bt+(b & 15u)*NT;
    double* hi = bt+(16+(b>>4))*NT;
    for(std::size_t k=0; k<NT; ++k){
      lo[k] += ti[k];
      hi[k] += ti[k];
    }
  };
  (one(std::integral_constant<std::size_t,L>()), ...);
}

// Every record of multTermsT at one column
template<std::size_t P, int M, std::size_t NT>
void multTermsTFixed(double* binPtr, const unsigned char* col,
                     const double* tPtr, arma::uword nInd){
  for(arma::uword i=0; i<nInd; ++i){
    multTermsTRecord<P,M,NT>(binPtr, col+i*P, tPtr+i*NT, i & 1u,
                             std::make_index_sequence<FixedLooks<P,M>::N>{});
  }
}

// The model of a list of terms, or -1 for one with no fixed loops
inline int fixedModel(const std::vector<FastTerm>& terms){
  if( (terms.size()==2) && (terms[0].kind==TermKind::additive) &&
      (terms[1].kind==TermKind::dominance) ){
    return 0;
  }
  if( (terms.size()==2) && (terms[0].kind==TermKind::maternal) &&
      (terms[1].kind==TermKind::paternal) ){
    return 1;
  }
  if( (terms.size()==3) && (terms[0].kind==TermKind::maternal) &&
      (terms[1].kind==TermKind::paternal) &&
      (terms[2].kind==TermKind::dominance) ){
    return 2;
  }
  return -1;
}

// Whether the products take the fixed loops: a diploid or tetraploid, with
// up to maxFixedTraits traits, fitting one of the three models
inline bool fixedTemplated(arma::uword ploidy, arma::uword nT,
                           const std::vector<FastTerm>& terms){
  return fixedPloidy(ploidy) && (nT>=1) && (nT<=maxFixedTraits) &&
    (fixedModel(terms)>=0);
}

// Calls f(std::integral_constant<int,M>()) for M equal to model, which
// fixedModel must have returned
template<class F>
void withFixedModel(int model, F&& f){
  if(model==0){
    f(std::integral_constant<int,0>());
  }else if(model==1){
    f(std::integral_constant<int,1>());
  }else{
    f(std::integral_constant<int,2>());
  }
}

// The general record loops of multTerms and multTermsT, for any ploidy and
// number of traits, in functions of their own, as multMcBlockGeneral and
// multMcTBlockGeneral are in MME.cpp, so that they do not share registers
// with the fixed ones in the parallel region.

// Records [indStart,indEnd) of multTerms
void multTermsGeneral(arma::mat& Tt, const arma::Mat<unsigned char>& planes,
                      arma::uword ploidy, const std::vector<FastTerm>& terms,
                      const arma::mat& tabs, arma::uword indStart,
                      arma::uword indEnd){
  const arma::uword nTerm = terms.size();
  const arma::uword nT = Tt.n_rows;
  const arma::uword nCol = planes.n_cols;
  const arma::uword termLen = 32*nT;
  for(arma::uword c=0; c<nCol; ++c){
    const unsigned char* col = planes.colptr(c);
    const double* tabCol = tabs.colptr(c);
    for(arma::uword i=indStart; i<indEnd; ++i){
      double* ti = Tt.colptr(i);
      const unsigned char* row = col+i*ploidy;
      for(arma::uword t=0; t<nTerm; ++t){
        const FastTerm& term = terms[t];
        const double* tab = tabCol+t*termLen;
        if(term.kind==TermKind::dominance){
          for(arma::uword h=0; h+1<ploidy; ++h){
            for(arma::uword h2=h+1; h2<ploidy; ++h2){
              addLookup(ti, tab, (unsigned int) (row[h]^row[h2]), nT);
            }
          }
        }else{
          for(arma::uword h=term.hapStart; h<term.hapEnd; ++h){
            addLookup(ti, tab, (unsigned int) row[h], nT);
          }
        }
      }
    }
  }
}

// Applies every term's centered code matrix at once,
//   Tt = (sum_k Mc_k*Ut_k')'
// with one row of Tt per trait and one column per record. offset holds, for
// each marker and term, the constant of the code plus the mean of the code,
// which centering takes off. This is multMcBlock with a table for each term,
// and as there the work is split over blocks of records, each of which owns
// its columns of Tt and sums the columns of the genotypes in order, so the
// result does not depend on the number of threads.
void multTerms(arma::mat& Tt, const arma::Mat<unsigned char>& planes,
               arma::uword ploidy, const std::vector<FastTerm>& terms,
               const arma::mat& offset, const std::vector<arma::mat>& Ut,
               arma::mat& tabs, int nThreads){
  arma::uword nTerm = terms.size();
  arma::uword nT = Tt.n_rows;
  arma::uword nInd = Tt.n_cols;
  arma::uword nCol = planes.n_cols;
  arma::uword nLoci = offset.n_rows;
  arma::uword termLen = 32*nT;
  arma::vec shift(nT, arma::fill::zeros);
  for(arma::uword t=0; t<nTerm; ++t){
    for(arma::uword j=0; j<nLoci; ++j){
      for(arma::uword k=0; k<nT; ++k){
        shift(k) += offset(j,t)*Ut[t](k,j);
      }
    }
    buildTables(tabs.memptr()+t*termLen, tabs.n_rows, Ut[t], terms[t].scale,
                nCol);
  }
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  // A diploid or tetraploid fitting one of the three models, with at most
  // maxFixedTraits traits, takes multTermsFixed
  const bool fixed = fixedTemplated(ploidy, nT, terms);
  const int model = fixedModel(terms);
  const double* tabsPtr = tabs.memptr();
  const arma::uword tabStride = tabs.n_rows;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword i=indStart; i<indEnd; ++i){
      double* ti = Tt.colptr(i);
      for(arma::uword k=0; k<nT; ++k){
        ti[k] = -shift(k);
      }
    }
    if(fixed){
      double* tPtr = Tt.memptr();
      withFixedPloidy(ploidy, [&](auto p){
        withFixedModel(model, [&](auto m){
          withFixedTraits(nT, [&](auto nt){
            multTermsFixed<decltype(p)::value, decltype(m)::value,
                           decltype(nt)::value>(tPtr, planes, tabsPtr,
                                                tabStride, indStart, indEnd);
          });
        });
      });
      continue;
    }
    multTermsGeneral(Tt, planes, ploidy, terms, tabs, indStart, indEnd);
  }
}

// Adds an individual's results into the bins of one byte
inline void addBins(double* bins, const double* ti, unsigned int b,
                    arma::uword nT){
  double* lo = bins+(b & 15u)*nT;
  double* hi = bins+(16+(b>>4))*nT;
  for(arma::uword k=0; k<nT; ++k){
    lo[k] += ti[k];
    hi[k] += ti[k];
  }
}

// Every record of multTermsT at one column, into each term's two sets of
// bins of 32*nT values
void multTermsTGeneral(double* binPtr, const unsigned char* col,
                       const arma::mat& Tt, arma::uword ploidy,
                       const std::vector<FastTerm>& terms){
  const arma::uword nTerm = terms.size();
  const arma::uword nT = Tt.n_rows;
  const arma::uword nInd = Tt.n_cols;
  const arma::uword setLen = 32*nT;
  for(arma::uword i=0; i<nInd; ++i){
    const double* ti = Tt.colptr(i);
    const unsigned char* row = col+i*ploidy;
    for(arma::uword t=0; t<nTerm; ++t){
      const FastTerm& term = terms[t];
      double* bt = binPtr+(2*t+(i & 1u))*setLen;
      if(term.kind==TermKind::dominance){
        for(arma::uword h=0; h+1<ploidy; ++h){
          for(arma::uword h2=h+1; h2<ploidy; ++h2){
            addBins(bt, ti, (unsigned int) (row[h]^row[h2]), nT);
          }
        }
      }else{
        for(arma::uword h=term.hapStart; h<term.hapEnd; ++h){
          addBins(bt, ti, (unsigned int) row[h], nT);
        }
      }
    }
  }
}

// Applies the transpose of every term's centered code matrix,
//   Wt_k = (Mc_k'*Tt')'
// one matrix per term, with one row per trait and one column per marker. As
// in multMcTBlock, each term's values are added into bins chosen by the bit
// pattern of each byte it reads, records alternating between two sets of
// bins that are merged in a fixed order, and a marker's sum is then the sum
// of the bins whose pattern has its bit set. Split over blocks of columns,
// each of which works out its own columns of every Wt_k, so nothing is
// accumulated across threads.
void multTermsT(std::vector<arma::mat>& Wt, const arma::Mat<unsigned char>& planes,
                arma::uword ploidy, const std::vector<FastTerm>& terms,
                const arma::mat& offset, const arma::mat& Tt, int nThreads){
  arma::uword nTerm = terms.size();
  arma::uword nT = Tt.n_rows;
  arma::uword nInd = Tt.n_cols;
  arma::uword nCol = planes.n_cols;
  arma::uword nLoci = offset.n_rows;
  arma::vec tSum(nT, arma::fill::zeros);
  for(arma::uword i=0; i<nInd; ++i){
    const double* ti = Tt.colptr(i);
    for(arma::uword k=0; k<nT; ++k){
      tSum(k) += ti[k];
    }
  }
  arma::uword nBlocks = countBlocks(nCol);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  // A diploid or tetraploid fitting one of the three models, with at most
  // maxFixedTraits traits, takes multTermsTFixed
  const bool fixed = fixedTemplated(ploidy, nT, terms);
  const int model = fixedModel(terms);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword colStart = blockStart(nCol, nBlocks, block);
    arma::uword colEnd = blockStart(nCol, nBlocks, block+1);
    // For each term two sets of 32 bins of nT values: in each, the low half
    // byte's bins, then the high's
    const arma::uword setLen = 32*nT;
    std::vector<double> bins(2*nTerm*setLen), merged(setLen);
    for(arma::uword c=colStart; c<colEnd; ++c){
      const unsigned char* col = planes.colptr(c);
      std::fill(bins.begin(), bins.end(), 0.0);
      if(fixed){
        const double* tPtr = Tt.memptr();
        double* binPtr = bins.data();
        withFixedPloidy(ploidy, [&](auto p){
          withFixedModel(model, [&](auto m){
            withFixedTraits(nT, [&](auto nt){
              multTermsTFixed<decltype(p)::value, decltype(m)::value,
                              decltype(nt)::value>(binPtr, col, tPtr, nInd);
            });
          });
        });
      }else{
        multTermsTGeneral(bins.data(), col, Tt, ploidy, terms);
      }
      for(arma::uword t=0; t<nTerm; ++t){
        const double* b0 = bins.data()+2*t*setLen;
        const double* b1 = b0+setLen;
        for(arma::uword e=0; e<setLen; ++e){
          merged[e] = b0[e]+b1[e];
        }
        double scale = terms[t].scale;
        for(arma::uword b=0; b<8; ++b){
          arma::uword j = 8*c+b;
          if(j>=nLoci){
            break;
          }
          const double* mb = merged.data()+16*(b/4)*nT;
          unsigned int bit = 1u<<(b%4);
          double* wj = Wt[t].colptr(j);
          for(arma::uword k=0; k<nT; ++k){
            double sum = 0.0;
            for(unsigned int q=1; q<16; ++q){
              if(q & bit){
                sum += mb[q*nT+k];
              }
            }
            wj[k] = scale*sum-offset(j,t)*tSum(k);
          }
        }
      }
    }
  }
}

// The dominance code summed over every marker, for each record, which the
// direct form of the models fits as a covariate for inbreeding depression.
// Bits past the last marker are zero on every haplotype, so they add nothing.
// Split over blocks of records, each summing its own records in order.
arma::vec dominanceSum(const arma::Mat<unsigned char>& planes,
                       arma::uword ploidy, double scale, int nThreads){
  arma::uword nInd = planes.n_rows/ploidy;
  arma::uword nCol = planes.n_cols;
  unsigned int popcount[256];
  for(unsigned int v=0; v<256; ++v){
    unsigned int n = 0;
    for(unsigned int b=v; b>0; b>>=1){
      n += b & 1u;
    }
    popcount[v] = n;
  }
  arma::vec out(nInd);
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
    for(arma::uword i=indStart; i<indEnd; ++i){
      // A count of pairs is an integer, so it is summed exactly
      arma::uword count = 0;
      for(arma::uword c=0; c<nCol; ++c){
        const unsigned char* row = planes.colptr(c)+i*ploidy;
        for(arma::uword h=0; h+1<ploidy; ++h){
          for(arma::uword h2=h+1; h2<ploidy; ++h2){
            count += popcount[row[h]^row[h2]];
          }
        }
      }
      out(i) = scale*double(count);
    }
  }
  return out;
}

// What one pass over the markers collects. Mmean holds the mean code of each
// term, one column per term. Dkl holds, for each marker, the diagonal entry
// of Mc_k'*S*Mc_l for every pair of terms, as a column of nTerm*nTerm values.
// genoCount holds the number of records with each dosage, and both11 the
// number carrying the allele on both haplotypes, for a diploid.
struct MarkerStats{
  arma::mat Mmean;
  arma::mat Dkl;
  arma::mat genoCount;
  arma::rowvec both11;
};

// The pass over the markers. Each marker is read once, its codes worked out
// for every term from the alleles, and everything accumulated into values
// for that marker alone, so the pass splits over blocks of markers with no
// memory shared between threads.
MarkerStats markerStats(const arma::Mat<unsigned char>& planes,
                        arma::uword ploidy, arma::uword nLoci,
                        const std::vector<FastTerm>& terms,
                        const arma::mat& Q, int nThreads){
  arma::uword nTerm = terms.size();
  arma::uword nInd = planes.n_rows/ploidy;
  arma::uword p = Q.n_cols;
  arma::rowvec Qsum = sum(Q,0); // Q'*1
  MarkerStats st;
  st.Mmean.set_size(nLoci,nTerm);
  st.Dkl.set_size(nTerm*nTerm,nLoci);
  st.genoCount.set_size(ploidy+1,nLoci);
  st.both11.set_size(nLoci);
  const double* QPtr = Q.memptr();
  double dn = double(nInd);
  arma::uword nBlocks = countBlocks(nLoci);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    std::vector<unsigned int> b(ploidy);
    arma::vec v(nTerm), s1(nTerm), mean(nTerm);
    arma::mat sKL(nTerm,nTerm), QtM(p,nTerm);
    std::vector<double> count(ploidy+1);
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      s1.zeros();
      sKL.zeros();
      QtM.zeros();
      std::fill(count.begin(), count.end(), 0.0);
      double nBoth = 0.0;
      for(arma::uword i=0; i<nInd; ++i){
        readBits(planes, ploidy, i, j, b.data());
        unsigned int x = 0;
        for(arma::uword h=0; h<ploidy; ++h){
          x += b[h];
        }
        count[x] += 1.0;
        if( (ploidy==2) && (b[0]==1u) && (b[1]==1u) ){
          nBoth += 1.0;
        }
        for(arma::uword t=0; t<nTerm; ++t){
          v(t) = termValue(terms[t], b.data(), ploidy);
        }
        for(arma::uword t=0; t<nTerm; ++t){
          s1(t) += v(t);
          for(arma::uword t2=t; t2<nTerm; ++t2){
            sKL(t,t2) += v(t)*v(t2);
          }
          for(arma::uword r=0; r<p; ++r){
            QtM(r,t) += QPtr[i+r*nInd]*v(t);
          }
        }
      }
      for(arma::uword t=0; t<nTerm; ++t){
        mean(t) = s1(t)/dn;
        st.Mmean(j,t) = mean(t);
        // Centering a column takes its mean times Q'*1 off Q'*M
        for(arma::uword r=0; r<p; ++r){
          QtM(r,t) -= mean(t)*Qsum(r);
        }
      }
      // Centering takes n times the product of the means off M_k'*M_l, and
      // absorbing the fixed effects takes off the product of the projections
      for(arma::uword t=0; t<nTerm; ++t){
        for(arma::uword t2=t; t2<nTerm; ++t2){
          double d = sKL(t,t2)-dn*mean(t)*mean(t2)-
            arma::dot(QtM.col(t), QtM.col(t2));
          if( (t==t2) && (d<0.0) ){
            // Rounding only, the true value cannot be negative
            d = 0.0;
          }
          st.Dkl(t*nTerm+t2,j) = d;
          st.Dkl(t2*nTerm+t,j) = d;
        }
      }
      for(arma::uword x=0; x<=ploidy; ++x){
        st.genoCount(x,j) = count[x];
      }
      st.both11(j) = nBoth;
    }
  }
  return st;
}

// A term's codes for a subset of records, as a dense matrix, for the REML
// fit of the variance components. With fewer markers than records in the
// subset these are the codes themselves. Otherwise it is the cross product of
// the codes, M*M', built a panel of markers at a time, which is all the
// record side of the REML fit uses. The model depends on a term only through
// M*M', so the variance components are the same either way and on the same
// scale, the variance of a marker effect.
arma::mat subsetCodes(const arma::Mat<unsigned char>& planes,
                      arma::uword ploidy, arma::uword nLoci,
                      const FastTerm& term, const arma::uvec& subset,
                      int nThreads){
  arma::uword nSub = subset.n_elem;
  bool useCodes = nLoci<nSub;
  arma::uword panelWidth = useCodes ? nLoci : 256;
  arma::mat out, K;
  if(useCodes){
    out.set_size(nSub,nLoci);
  }else{
    K.zeros(nSub,nSub);
  }
  for(arma::uword start=0; start<nLoci; start+=panelWidth){
    arma::uword stop = std::min(start+panelWidth, nLoci);
    arma::mat panel(nSub,stop-start);
    arma::uword nBlocks = countBlocks(stop-start);
    int panelThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(panelThreads)){
      panelThreads = static_cast<int>(nBlocks);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(panelThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword jStart = start+blockStart(stop-start, nBlocks, block);
      arma::uword jEnd = start+blockStart(stop-start, nBlocks, block+1);
      std::vector<unsigned int> b(ploidy);
      for(arma::uword j=jStart; j<jEnd; ++j){
        double* col = panel.colptr(j-start);
        for(arma::uword i=0; i<nSub; ++i){
          readBits(planes, ploidy, subset(i), j, b.data());
          col[i] = termValue(term, b.data(), ploidy);
        }
      }
    }
    if(useCodes){
      out = panel;
    }else{
      K += panel*panel.t();
    }
  }
  if(useCodes){
    return out;
  }
  return K;
}

// A thin factor F of a term's cross product K = M*M' on the subset, with
// F*F' = K and as many columns as the numerical rank of K. The rank matters
// twice. The EM step of the REML fit divides each update by a bound for it,
// and the number of records, the only bound otherwise to hand, is far too
// loose for a kernel of low rank, as from many copies of a few genotypes.
// And a kernel of low rank is cheaper to fit on the coefficient side, from
// its factor, than on the record side from K.
//
// The factor comes from a pivoted Cholesky factorization of a copy of K,
// P'*K*P = L*L' with L lower trapezoidal, so F = P*L. It costs about a third
// of a cube of the number of records, against the four or more of a
// decomposition with eigenvectors, which a factor used to come from. The
// tolerance is that of rankCutoff, taken against the largest diagonal
// element, which is no larger than the largest eigenvalue, so the rank found
// is at least the number of eigenvalues that rankCutoff keeps. Returns false
// if the factorization is not available or fails, and then F is empty.
bool gramFactor(const arma::mat& K, arma::mat& F){
  arma::uword n = K.n_rows;
  F.reset();
  if(n==0){
    return false;
  }
#if ARMA_VERSION_MAJOR >= 10
  arma::mat A = K;
  char uplo = 'L';
  arma::blas_int bn = arma::blas_int(n);
  arma::blas_int rank = 0;
  arma::blas_int info = 0;
  double tol = rankCutoff(A.diag().max(), n);
  std::vector<arma::blas_int> piv(n);
  std::vector<double> work(2*n);
  arma::lapack::pstrf(&uplo, &bn, A.memptr(), &bn, piv.data(), &rank, &tol,
                      work.data(), &info);
  // info is 1 when K is rank deficient, which is the case that matters
  if( (info!=0) && (info!=1) ){
    return false;
  }
  arma::uword r = arma::uword(std::max<arma::blas_int>(rank, 0));
  // The factor is the first r columns of the lower triangle. Above the
  // diagonal A still holds K, and below the first r columns, the part of K
  // left unfactored, so both are left out.
  F.zeros(n, r);
  for(arma::uword j=0; j<r; ++j){
    for(arma::uword i=j; i<n; ++i){
      F(arma::uword(piv[i]-1), j) = A(i, j);
    }
  }
  return true;
#else
  return false;
#endif
}

// The result of a fast multikernel fit. U holds each term's marker effects,
// one column per trait, and beta the fixed effects for the uncentered codes,
// as the direct form returns them.
struct FastMKFit{
  std::vector<arma::mat> U;
  arma::mat beta;
  std::vector<arma::mat> Vu;
  arma::mat Ve;
  MarkerStats stats;
  arma::uword rank = 0; // eigenvalues the preconditioner deflated
};

// Sets fit.Vu and fit.Ve. VuIn and VeIn are used when estVarComp is false,
// as nT by nT matrices with one for each term. Otherwise they are estimated
// by REML on the records in subset, which is one based. As in the fast mode
// of RRBLUP, the decomposition REML needs costs the square of the number of
// records in memory and their cube in time, so it is done on a subset, and
// the marker effects are then solved for on every record.
void fastMKVarComp(FastMKFit& fit, const arma::mat& Y, const arma::mat& X,
                   const arma::Mat<unsigned char>& planes,
                   arma::uword ploidy, arma::uword nLoci,
                   const std::vector<FastTerm>& terms,
                   const arma::field<arma::mat>& VuIn,
                   const arma::mat& VeIn, bool estVarComp,
                   arma::uvec subset, int nThreads,
                   const std::string& name){
  const arma::uword nT = Y.n_cols;
  const arma::uword nTerm = terms.size();
  fit.Vu.resize(nTerm);
  if(estVarComp){
    subset -= 1; // R to C++
    arma::field<arma::mat> Mlist(nTerm);
    for(arma::uword t=0; t<nTerm; ++t){
      Mlist(t) = subsetCodes(planes, ploidy, nLoci, terms[t], subset,
                             nThreads);
    }
    arma::mat ySub = Y.rows(subset);
    arma::mat XSub = X.rows(subset);
    // With at least as many markers as records, subsetCodes gives each
    // term's cross product. Each is factored to find its rank, and the side
    // is chosen as the REML fit would choose it for factors with that many
    // columns. Usually the ranks are close to the number of records and the
    // cross products go to the record side directly. When they are low, as
    // for many copies of a few genotypes, the factors are fitted instead,
    // which the REML fit then does on the far smaller coefficient side.
    Rcpp::List varComp;
    const arma::uword nSub = subset.n_elem;
    if(nLoci<nSub){
      varComp = solveRRBLUPMK(ySub, XSub, Mlist, fastMKRemlIter);
    }else{
      std::vector<arma::mat> F(nTerm);
      arma::vec rank(nTerm);
      bool factored = true;
      arma::uword nCoef = XSub.n_cols;
      for(arma::uword t=0; t<nTerm; ++t){
        if(gramFactor(Mlist(t), F[t])){
          rank(t) = double(std::max<arma::uword>(F[t].n_cols, 1));
          nCoef += F[t].n_cols;
        }else{
          factored = false;
          rank(t) = double(nSub);
          nCoef += nSub;
        }
      }
      bool coefficientSide = factored &&
        ((nT==1) ? useCoefficientSide(nSub, nCoef, nTerm) :
           useCoefficientSide(nSub*nT, nCoef*nT,
                              (nTerm+1)*nT*(nT+1)/2));
      if(coefficientSide){
        for(arma::uword t=0; t<nTerm; ++t){
          Mlist(t) = std::move(F[t]);
        }
        varComp = solveRRBLUPMK(ySub, XSub, Mlist, fastMKRemlIter);
      }else{
        F.clear();
        varComp = remlMultiKernelGram(ySub, XSub, Mlist, rank,
                                      fastMKRemlIter, name.c_str());
      }
    }
    if(nT==1){
      arma::mat Vu = varComp["Vu"];
      for(arma::uword t=0; t<nTerm; ++t){
        fit.Vu[t].set_size(1,1);
        fit.Vu[t](0,0) = Vu(t);
      }
    }else{
      Rcpp::List Vu = varComp["Vu"];
      for(arma::uword t=0; t<nTerm; ++t){
        fit.Vu[t] = Rcpp::as<arma::mat>(Vu[t]);
      }
    }
    fit.Ve = Rcpp::as<arma::mat>(varComp["Ve"]);
  }else{
    if(VuIn.n_elem!=nTerm){
      Rcpp::stop("%s: Vu must have one entry per random term", name.c_str());
    }
    for(arma::uword t=0; t<nTerm; ++t){
      fit.Vu[t] = 0.5*(VuIn(t)+VuIn(t).t());
    }
    fit.Ve = 0.5*(VeIn+VeIn.t());
  }
  for(arma::uword t=0; t<nTerm; ++t){
    if( (fit.Vu[t].n_rows!=nT) || (fit.Vu[t].n_cols!=nT) ){
      Rcpp::stop("%s: Vu and Ve must be square with one row per trait",
                 name.c_str());
    }
  }
  if( (fit.Ve.n_rows!=nT) || (fit.Ve.n_cols!=nT) ){
    Rcpp::stop("%s: Vu and Ve must be square with one row per trait",
               name.c_str());
  }
  // Each Vu is checked on its own scale, as in the fast mode of RRBLUP,
  // because after whitening its eigenvalues are measured against those of
  // Ve^-1. A negative eigenvalue at rounding level is set to zero, and the
  // Vu returned is the one fitted.
  std::string psdMsg = name+": Vu must be positive semidefinite";
  for(arma::uword t=0; t<nTerm; ++t){
    fit.Vu[t] = psdKernel(fit.Vu[t], psdMsg);
  }
}

// The model on the scale whitened by Ve = Le*Le'. Each term's covariance
// there is F*F', and its whitened effects for term t are the rows
// off[t] to off[t+1]-1 of the stacked effects. Finv maps a term's effects
// back to its whitened effects, and Fdir holds the trait directions the term
// has variance in.
struct Whitened{
  arma::mat Le;
  arma::mat Linv;
  std::vector<arma::mat> F;
  std::vector<arma::mat> Finv;
  std::vector<arma::mat> Fdir;
  std::vector<arma::uword> off;
};

Whitened whitenTerms(const std::vector<arma::mat>& Vu, const arma::mat& Ve,
                     const std::string& name){
  const arma::uword nTerm = Vu.size();
  Whitened w;
  if(!arma::chol(w.Le, Ve, "lower")){
    Rcpp::stop("%s: Ve must be positive definite", name.c_str());
  }
  if(!arma::inv(w.Linv, arma::trimatl(w.Le))){
    Rcpp::stop("%s: Ve must be positive definite", name.c_str());
  }
  // These are eigenvalues of a covariance matrix among traits, so rankCutoff
  // does not apply, and a positive eigenvalue is kept however small it is
  // next to the others, for the reason given in callFastRRBLUP_MV. A tiny one
  // only adds a direction that barely moves.
  w.F.resize(nTerm);
  w.Finv.resize(nTerm);
  w.Fdir.resize(nTerm);
  w.off.assign(nTerm+1, 0);
  for(arma::uword t=0; t<nTerm; ++t){
    arma::mat A = w.Linv*Vu[t]*w.Linv.t();
    A = 0.5*(A+A.t());
    arma::vec lam;
    arma::mat Qv;
    if(!arma::eig_sym(lam, Qv, A)){
      Rcpp::stop("%s: decomposition of Vu failed", name.c_str());
    }
    arma::uvec keep = find(lam>0.0);
    w.F[t] = Qv.cols(keep);
    w.F[t].each_row() %= sqrt(lam(keep)).t();
    w.Finv[t] = Qv.cols(keep).t();
    w.Finv[t].each_col() /= sqrt(lam(keep));
    w.Fdir[t] = Qv.cols(keep);
    w.off[t+1] = w.off[t]+keep.n_elem;
  }
  return w;
}

// The preconditioner is block Jacobi, with one block for each marker over
// every term and whitened trait. The terms read at one marker are often
// correlated, as the additive and dominance codes are wherever the allele
// frequency is not one half, so the blocks hold their cross products as
// well as the diagonal. Each block is the identity plus a positive
// semidefinite matrix G, so it always has an inverse, but finding it takes
// care, because the entries of a block can differ in scale by more than
// the precision of a double. They do when one trait's error variance is
// far below another's, which whitening turns into variance ratios many
// orders of magnitude apart, and when terms coincide at a marker, as the
// maternal and paternal codes of inbred lines do, the identity is lost to
// rounding next to G and the block is singular in floating point.
//
// The block is therefore first scaled to unit diagonal, Ds*B*Ds for
// Ds = diag(1/sqrt(1+diag(G))), and its eigendecomposition taken. The
// scaled block's entries are all of order one, so every direction the
// records see has an eigenvalue found to rounding level, where the
// unscaled block's decomposition would find small ones only to rounding of
// its largest entry, mixing them up with each other and with the unseen
// directions. An eigenvalue at rounding level is not resolved by the
// block: it may belong to a direction the records do not see, which the
// deflation removes, or to one they see only through differences
// below the block's precision, as terms whose covariance matrices are
// nearly but not exactly aligned have. Its inverse is taken as one, the
// scaled block's diagonal, rather than the reciprocal of rounding error or
// zero: the preconditioner then stays positive definite, so every
// direction is still searched in, and bounded. The inverse is kept and
// applied in factored form, Ds*V*(dInv%(V'*(Ds*r))),
// rather than multiplied out, so that its scales are never added together.
// Should the decomposition fail, the inverse of the diagonal is used, which
// is still a positive definite preconditioner.
struct BlockPrecond{
  arma::cube V;   // eigenvectors of each marker's scaled block
  arma::mat dInv; // inverses of their eigenvalues, one column per marker
  arma::mat ds;   // the diagonal of Ds, one column per marker
};

BlockPrecond buildPrecond(const MarkerStats& stats, const Whitened& w,
                          arma::uword nLoci, int nThreads){
  const arma::uword nTerm = w.F.size();
  const std::vector<arma::uword>& off = w.off;
  const arma::uword nR = off[nTerm];
  BlockPrecond P;
  P.V.set_size(nR,nR,nLoci);
  P.dInv.set_size(nR,nLoci);
  P.ds.set_size(nR,nLoci);
  if(nR==0){
    return P;
  }
  arma::uword nBlocks = countBlocks(nLoci);
  int precondThreads = nThreads;
  if(nBlocks < static_cast<arma::uword>(precondThreads)){
    precondThreads = static_cast<int>(nBlocks);
  }
  // F_k'*F_l for every pair of terms, the same for every marker
  std::vector<arma::mat> FF(nTerm*nTerm);
  for(arma::uword t=0; t<nTerm; ++t){
    for(arma::uword t2=0; t2<nTerm; ++t2){
      FF[t*nTerm+t2] = w.F[t].t()*w.F[t2];
    }
  }
  // An exception cannot leave the parallel region, see unseenDirections
  std::vector<int> failed(nBlocks, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(precondThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    try{
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      arma::mat B(nR,nR), V;
      arma::vec lambda;
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        B.zeros();
        for(arma::uword t=0; t<nTerm; ++t){
          for(arma::uword t2=0; t2<nTerm; ++t2){
            if( (off[t+1]==off[t]) || (off[t2+1]==off[t2]) ){
              continue;
            }
            B.submat(off[t], off[t2], off[t+1]-1, off[t2+1]-1) +=
              stats.Dkl(t*nTerm+t2,j)*FF[t*nTerm+t2];
          }
        }
        B = 0.5*(B+B.t());
        arma::vec ds = 1.0/sqrt(1.0+arma::clamp(B.diag(), 0.0,
                                                arma::datum::inf));
        P.ds.col(j) = ds;
        B.each_col() %= ds;
        B.each_row() %= ds.t();
        B.diag() += square(ds);
        if(arma::eig_sym(lambda, V, B)){
          double cut = 64.0*arma::datum::eps*std::max(lambda.max(), 0.0);
          P.V.slice(j) = V;
          arma::vec dInv(nR, arma::fill::ones);
          for(arma::uword c=0; c<nR; ++c){
            if(lambda(c)>cut){
              dInv(c) = 1.0/lambda(c);
            }
          }
          P.dInv.col(j) = dInv;
        }else{
          P.V.slice(j).eye();
          P.dInv.col(j).ones();
        }
      }
    }catch(...){
      failed[block] = 1;
    }
  }
  if(std::find(failed.begin(), failed.end(), 1)!=failed.end()){
    Rcpp::stop("fast mode: building the preconditioner failed");
  }
  return P;
}

// Z = P^-1*R, each marker's block applied to its own column, so the work
// splits over blocks of markers
void applyPrecond(arma::mat& Z, const arma::mat& R, const BlockPrecond& P,
                  int nThreads){
  const arma::uword nR = R.n_rows;
  const arma::uword nLoci = R.n_cols;
  arma::uword nBlocks = countBlocks(nLoci);
  int precondThreads = nThreads;
  if(nBlocks < static_cast<arma::uword>(precondThreads)){
    precondThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(precondThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    std::vector<double> w(nR);
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const double* Vj = P.V.slice_memptr(j);
      const double* dj = P.dInv.colptr(j);
      const double* sj = P.ds.colptr(j);
      const double* rj = R.colptr(j);
      double* zj = Z.colptr(j);
      // w = dInv%(V'*(Ds*r)), then z = Ds*(V*w)
      for(arma::uword c=0; c<nR; ++c){
        const double* vc = Vj+c*nR;
        double sum = 0.0;
        for(arma::uword a=0; a<nR; ++a){
          sum += vc[a]*(sj[a]*rj[a]);
        }
        w[c] = dj[c]*sum;
      }
      for(arma::uword a=0; a<nR; ++a){
        double sum = 0.0;
        for(arma::uword c=0; c<nR; ++c){
          sum += Vj[a+c*nR]*w[c];
        }
        zj[a] = sj[a]*sum;
      }
    }
  }
}

// A deflation of the block preconditioner, as the fast mode of RRBLUP
// deflates its Jacobi preconditioner; see deflationBasis in MME.cpp.
//
// The coefficient matrix is A = I+G'*G, where G takes the stacked whitened
// effects to the whitened fitted values of every record and trait, and the
// block preconditioner's inverse is L*L', with L_j = Ds_j*V_j*dInv_j^1/2 for
// marker j. The iteration then works, in effect, on L'*A*L, which has a mean
// eigenvalue near one but, as for one trait, a few hundred to a few
// thousand eigenvalues far above it from the relationships between
// individuals, and those set the number of iterations. If W holds
// orthonormal approximations to their eigenvectors and theta their values,
//   P^-1 = L*(I+W*diag(sigma/theta-1)*W')*L'
// moves them down to about sigma and leaves the rest alone, and is positive
// definite for any orthonormal W, so a rough W costs iterations, never the
// solution.
//
// W is estimated from the rows of G*L for nSub records evenly spaced through
// the population, with their codes centered on the subset, and every trait,
// so nSub*nT rows. n/nSub*(Gs*L)'*(Gs*L) approximates (G*L)'*(G*L), and the
// mean of the diagonal of L'*L stands in for L'*L on every eigenvalue. The
// leading eigenvectors are (Gs*L)'*U*diag(1/sqrt(g)) for the eigenpairs
// (g,U) of the cross product of the rows, (Gs*L)*(Gs*L)', which is formed a
// panel of markers at a time. As for one trait, the fixed effects other
// than the mean are left out of the approximation.
//
// The records are chosen so that nSub*nT is about nPrecond, which keeps the
// decomposition the size the fast mode of RRBLUP uses for one trait. W has a
// row for every whitened direction of every marker, nR of them, so at most
// maxRank/nR eigenvalues are deflated, which holds W to the memory the fast
// mode of RRBLUP gives it. Eigenvalues above two are deflated and sigma is
// the largest one left, never less than one, as there.
//
// The directions the records cannot see are left to deflateUnseen: they are
// in the null space of G*L once scaled, so W, built from the rows of Gs*L,
// is close to orthogonal to them, and the projection that follows the
// preconditioner removes what is left.
struct MKDeflation{
  arma::cube L;     // L_j for each marker
  arma::mat W;      // row j*nR+c for direction c of marker j
  arma::vec shrink; // sigma/theta-1
};

MKDeflation mkDeflation(const arma::Mat<unsigned char>& planes,
                        arma::uword ploidy, const std::vector<FastTerm>& terms,
                        const Whitened& w, const BlockPrecond& P,
                        arma::uword nT, arma::uword nInd, arma::uword nLoci,
                        arma::uword nPrecond, arma::uword maxRank,
                        int nThreads){
  MKDeflation D;
  const arma::uword nTerm = terms.size();
  const std::vector<arma::uword>& off = w.off;
  const arma::uword nR = off[nTerm];
  if((nR==0) || (nPrecond==0) || (maxRank==0) || (nT==0)){
    return D;
  }
  arma::uword nSub = std::max<arma::uword>(nPrecond/nT, 1);
  if(nSub>nInd){
    nSub = nInd;
  }
  const arma::uword nRow = nSub*nT;
  maxRank /= nR;
  arma::uword nDecomp = std::min(nRow, nR*nLoci);
  if(maxRank>=nDecomp){
    maxRank = (nDecomp>0) ? nDecomp-1 : 0;
  }
  if(maxRank==0){
    return D;
  }

  // The factor of each marker's block inverse
  D.L.set_size(nR, nR, nLoci);
  {
    arma::uword nBlocks = countBlocks(nLoci);
    int lThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(lThreads)){
      lThreads = static_cast<int>(nBlocks);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(lThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        const double* V = P.V.slice_memptr(j);
        const double* dInv = P.dInv.colptr(j);
        const double* ds = P.ds.colptr(j);
        double* Lj = D.L.slice_memptr(j);
        for(arma::uword c=0; c<nR; ++c){
          double s = std::sqrt(dInv[c]);
          for(arma::uword a=0; a<nR; ++a){
            Lj[a+c*nR] = ds[a]*V[a+c*nR]*s;
          }
        }
      }
    }
  }
  // The mean of the diagonal of L'*L
  double ridge = arma::accu(arma::square(D.L))/double(nR*nLoci);

  arma::uvec sub = deflationSubset(nInd, nSub);
  const arma::uword panelWidth = std::max<arma::uword>(
    16, 8388608/std::max<arma::uword>(nRow*nR, 1));
  // One panel of Gs*L, for markers [start,stop), with row k*nSub+i for
  // record i and whitened trait k and column (j-start)*nR+c
  std::vector<unsigned int> bits(ploidy);
  arma::mat codes(nSub, nTerm), H;
  auto fillPanel = [&](arma::mat& panel, arma::uword start, arma::uword stop){
    panel.zeros(nRow, (stop-start)*nR);
    for(arma::uword j=start; j<stop; ++j){
      for(arma::uword i=0; i<nSub; ++i){
        readBits(planes, ploidy, sub(i), j, bits.data());
        for(arma::uword t=0; t<nTerm; ++t){
          codes(i,t) = termValue(terms[t], bits.data(), ploidy);
        }
      }
      codes.each_row() -= arma::mean(codes, 0);
      const arma::mat Lj(const_cast<double*>(D.L.slice_memptr(j)), nR, nR,
                         false, true);
      for(arma::uword t=0; t<nTerm; ++t){
        if(off[t+1]==off[t]){
          continue;
        }
        // Term t's codes reach whitened trait k through F_t, and the scaled
        // direction c through the rows of L_j for the term
        H = w.F[t]*Lj.rows(off[t], off[t+1]-1);
        const double* ct = codes.colptr(t);
        for(arma::uword c=0; c<nR; ++c){
          double* pc = panel.colptr((j-start)*nR+c);
          for(arma::uword k=0; k<nT; ++k){
            double h = H(k,c);
            double* pk = pc+k*nSub;
            for(arma::uword i=0; i<nSub; ++i){
              pk[i] += ct[i]*h;
            }
          }
        }
      }
    }
  };

  arma::mat K(nRow, nRow, arma::fill::zeros), panel;
  for(arma::uword start=0; start<nLoci; start+=panelWidth){
    arma::uword stop = std::min(start+panelWidth, nLoci);
    fillPanel(panel, start, stop);
    K += panel*panel.t();
  }
  arma::vec g;
  arma::mat U;
  if(!arma::eig_sym(g, U, K)){
    // The block preconditioner alone still converges
    return MKDeflation();
  }
  K.reset();
  // eig_sym returns the eigenvalues in ascending order
  g = arma::flipud(g);
  U = arma::fliplr(U);
  double scale = double(nInd)/double(nSub);
  arma::vec theta = scale*g+ridge;
  arma::uword k = 0;
  while((k<maxRank) && (theta(k)>2.0) && (g(k)>0.0)){
    ++k;
  }
  if(k==0){
    return MKDeflation();
  }
  double sigma = std::max(theta(k), 1.0);
  arma::mat Uk = U.head_cols(k);
  U.reset();
  Uk.each_row() /= arma::sqrt(g.head(k)).t();
  D.W.set_size(nR*nLoci, k);
  for(arma::uword start=0; start<nLoci; start+=panelWidth){
    arma::uword stop = std::min(start+panelWidth, nLoci);
    fillPanel(panel, start, stop);
    D.W.rows(start*nR, stop*nR-1) = panel.t()*Uk;
  }
  D.shrink = sigma/theta.head(k)-1.0;
  return D;
}

// Z = P^-1*R with the deflation, L*(I+W*diag(shrink)*W')*L'*R, each marker's
// factor applied to its own column, so the work splits over blocks of
// markers, and the deflation applied as for one trait, so neither depends on
// the number of threads
void applyDeflatedPrecond(arma::mat& Z, const arma::mat& R,
                          const MKDeflation& D, int nThreads){
  const arma::uword nR = R.n_rows;
  const arma::uword nLoci = R.n_cols;
  arma::vec v(nR*nLoci);
  arma::uword nBlocks = countBlocks(nLoci);
  int pThreads = nThreads;
  if(nBlocks < static_cast<arma::uword>(pThreads)){
    pThreads = static_cast<int>(nBlocks);
  }
  double* vPtr = v.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(pThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const double* Lj = D.L.slice_memptr(j);
      const double* rj = R.colptr(j);
      double* vj = vPtr+j*nR;
      for(arma::uword c=0; c<nR; ++c){
        const double* lc = Lj+c*nR;
        double sum = 0.0;
        for(arma::uword a=0; a<nR; ++a){
          sum += lc[a]*rj[a];
        }
        vj[c] = sum;
      }
    }
  }
  applyDeflation(v, D.W, D.shrink, nThreads);
  Z.set_size(nR, nLoci);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(pThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const double* Lj = D.L.slice_memptr(j);
      const double* vj = vPtr+j*nR;
      double* zj = Z.colptr(j);
      for(arma::uword a=0; a<nR; ++a){
        double sum = 0.0;
        for(arma::uword c=0; c<nR; ++c){
          sum += Lj[a+c*nR]*vj[c];
        }
        zj[a] = sum;
      }
    }
  }
}

// Directions the records cannot see. When some combination of the terms'
// codes at a marker is zero once centered and absorbed, as the maternal
// and paternal codes are for inbred lines, whose two haplotypes are the
// same, the effects along that combination do not change any prediction,
// and their BLUP is exactly zero: only the prior holds them. With a large
// variance ratio the data dwarf the prior, so rounding in the residual is
// as large as the prior's pull, and the iteration would let those effects
// drift. Their sum would stay right while the split between the terms went
// wrong. They are therefore removed exactly: every residual, search
// direction and product is kept orthogonal to them.
//
// The combinations are the null space of the terms' cross products at the
// marker, Dkl, among the terms with any variance. Those are all on the
// scale of the number of records, so their null space is found reliably,
// unlike that of the preconditioner's block, whose scales can differ by
// more than the precision of a double. For a null vector c, the effects
// F_k*a_k = c_k*w are unseen, for any w over the traits that every term in
// the combination can carry, and in whitened effects they are
// a_k = F_k^+*c_k*w. When every such term has variance in every trait
// direction, w is any vector and the traits' own directions serve as its
// basis. When some do not, as when a covariance matrix is singular, w is
// confined to the directions they all have variance in, the intersection
// of the column spaces of their F_k, found as the null space of the
// stacked projections onto the complements of those spaces. With more than
// one null vector the unseen effects are F_k*a_k = sum_q c_qk*w_q, and the
// weights are found together, because a combination of null vectors can be
// carried by every term when no single one of them is: the condition on
// vec(W) for each term k is (c_k'(x)P_k)*vec(W) = 0, for P_k the
// projection onto the trait directions it has no variance in.
//
// Returns an orthonormal basis of each marker's unseen directions in the
// stacked whitened effects, empty where there are none, and sets anyNull
// when any marker has one.
std::vector<arma::mat> unseenDirections(const MarkerStats& stats,
                                        const Whitened& w, arma::uword nT,
                                        arma::uword nInd, arma::uword nLoci,
                                        int nThreads, bool& anyNull){
  const arma::uword nTerm = w.F.size();
  const std::vector<arma::uword>& off = w.off;
  const arma::uword nR = off[nTerm];
  std::vector<arma::mat> nullBasis(nLoci);
  std::vector<arma::uword> active;
  for(arma::uword t=0; t<nTerm; ++t){
    if(off[t+1]>off[t]){
      active.push_back(t);
    }
  }
  arma::uword nA = active.size();
  arma::uword nBlocks = countBlocks(nLoci);
  int nullThreads = nThreads;
  if(nBlocks < static_cast<arma::uword>(nullThreads)){
    nullThreads = static_cast<int>(nBlocks);
  }
  std::vector<int> found(nBlocks, 0);
  // Anything else that throws inside the parallel region, such as an
  // allocation that fails, would end the R session, so it is caught and
  // reported after the region instead
  std::vector<int> failed(nBlocks, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nullThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    try{
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      arma::mat D(nA,nA), C, N;
      arma::vec dl;
      for(arma::uword j=lociStart; (j<lociEnd)&&(nA>0); ++j){
        for(arma::uword a=0; a<nA; ++a){
          for(arma::uword b=0; b<nA; ++b){
            D(a,b) = stats.Dkl(active[a]*nTerm+active[b],j);
          }
        }
        if(!arma::eig_sym(dl, C, D)){
          continue;
        }
        arma::uvec nul = find(dl<=rankCutoff(dl.max(), nInd));
        if(nul.n_elem==0){
          continue;
        }
        std::vector<arma::vec> cols;
        arma::mat Cn = C.cols(nul);
        arma::uword nq = nul.n_elem;
        double cMax = arma::abs(Cn).max();
        std::vector<arma::uword> involved;
        bool allFull = true;
        for(arma::uword a=0; a<nA; ++a){
          if(arma::abs(Cn.row(a)).max()>1.0e-8*cMax){
            involved.push_back(a);
            if(off[active[a]+1]-off[active[a]]!=nT){
              allFull = false;
            }
          }
        }
        // A basis for the weights vec(W), the trait values of every null
        // vector stacked
        arma::mat Wb;
        if(allFull){
          Wb.eye(nT*nq,nT*nq);
        }else{
          arma::mat E(0,nT*nq);
          for(arma::uword a : involved){
            arma::uword t = active[a];
            if(off[t+1]-off[t]==nT){
              continue;
            }
            arma::mat P = arma::eye(nT,nT)-w.Fdir[t]*w.Fdir[t].t();
            E = join_cols(E, arma::kron(Cn.row(a), P));
          }
          // The form that reports failure rather than throwing, because
          // an exception cannot leave a parallel region. The marker is
          // then left without a null basis, as when eig_sym fails.
          if(!arma::null(Wb, E)){
            continue;
          }
        }
        for(arma::uword i=0; i<Wb.n_cols; ++i){
          arma::mat W = arma::reshape(Wb.col(i), nT, nq);
          arma::vec col(nR, arma::fill::zeros);
          for(arma::uword a : involved){
            arma::uword t = active[a];
            col.subvec(off[t], off[t+1]-1) = w.Finv[t]*(W*Cn.row(a).t());
          }
          if(arma::norm(col)>0.0){
            cols.push_back(col);
          }
        }
        if(!cols.empty()){
          N.set_size(nR, cols.size());
          for(arma::uword i=0; i<cols.size(); ++i){
            N.col(i) = cols[i];
          }
          // The columns can differ in length by many orders of magnitude, as
          // F_k^-1 does, and an orthonormal basis taken through the singular
          // value decomposition would find the short ones' directions only to
          // the ratio of the lengths times the precision, mixing in part of
          // the directions the records do see. Each column is normalized
          // first, and the basis then taken by QR.
          for(arma::uword c=0; c<N.n_cols; ++c){
            N.col(c) /= arma::norm(N.col(c));
          }
          arma::mat Qn, Rn;
          if(arma::qr_econ(Qn, Rn, N)){
            nullBasis[j] = Qn;
            found[block] = 1;
          }
        }
      }
    }catch(...){
      failed[block] = 1;
    }
  }
  if(std::find(failed.begin(), failed.end(), 1)!=failed.end()){
    Rcpp::stop("fast mode: finding the unseen effects of a marker failed");
  }
  anyNull = std::find(found.begin(), found.end(), 1)!=found.end();
  return nullBasis;
}

// Removes the unseen directions from each marker's column of R
void deflateUnseen(arma::mat& R, const std::vector<arma::mat>& nullBasis,
                   bool anyNull){
  if(!anyNull){
    return;
  }
  for(arma::uword j=0; j<nullBasis.size(); ++j){
    const arma::mat& N = nullBasis[j];
    if(N.n_cols>0){
      R.col(j) -= N*(N.t()*R.col(j));
    }
  }
}

// Takes each level's first response off every record of the level, giving
// Yc, and returns the responses taken off, one row per level. A response
// with a large mean keeps its variation in its low order bits, which any
// product with it rounds away: whitening by Le^-1, and absorbing the fixed
// effects through Q, whose entries are not exact either. That subtraction is
// exact wherever the level's responses share a large offset, since two
// doubles within a factor of two of each other subtract without rounding,
// and it changes nothing else, because the level's effect absorbs a
// constant. Taking one record's response off every level instead would fail
// when the levels' means differ by more than the precision of a double, as
// the largest would round the others' variation away.
arma::mat levelOffsets(const arma::mat& Y, const arma::uvec& x,
                       arma::mat& Yc){
  const arma::uword nInd = Y.n_rows;
  arma::uword nLevel = (nInd>0) ? x.max() : 0;
  arma::uvec firstOfLevel(nLevel);
  std::vector<bool> seen(nLevel, false);
  for(arma::uword i=0; i<nInd; ++i){
    arma::uword l = x(i)-1;
    if(!seen[l]){
      seen[l] = true;
      firstOfLevel(l) = i;
    }
  }
  arma::mat levelOffset = Y.rows(firstOfLevel);
  Yc = Y;
  for(arma::uword i=0; i<nInd; ++i){
    Yc.row(i) -= levelOffset.row(x(i)-1);
  }
  return levelOffset;
}

// How restartedPCG stopped. stalledAt is the residual relative to the right
// hand sides when the iteration stalled, and roundingBound the last
// correction relative to the effects when rounding limits their accuracy;
// each is zero otherwise.
struct PcgStatus{
  bool converged;
  double stalledAt;
  double roundingBound;
};

// Preconditioned conjugate gradient on the stacked whitened effects, solving
// for At from the right hand side Bt. applyCoef(out, in) is one product with
// the coefficient matrix and precondition(Z, R) one with the inverse of the
// preconditioner, and lamMax is the largest diagonal entry of the
// coefficient matrix. The terms and traits share one recurrence, because the
// terms are coupled through the records.
//
// Whitening gives every error unit variance, but the right hand sides of
// different directions can still differ in scale by many orders of
// magnitude, as they do when one trait's error variance is close to zero,
// because F_k carries the variance ratio of each direction. Two things
// follow. A test on the residual as a whole would stop once the largest
// direction had converged, so each row of the residual, one whitened
// direction of one term, is tested against its own initial value, as if
// it were fitted on its own, down to rounding level of the largest. And the
// step lengths, which come from sums over every direction, are then
// decided by the largest, so the smaller directions' part of the
// recurrence is lost to rounding and its residual stops describing them.
// When the recurrence claims convergence, the true residual b-A*x is
// therefore worked out, and if any direction has not converged the
// iteration is restarted from it. By then the largest directions have
// converged, so the restarted sums are no longer dominated by them. Every
// product with the coefficient matrix counts toward maxIter.
//
// A residual that is small next to the right hand sides does not make the
// effects accurate, though. Directions that only the prior holds, such as
// effects that cancel across markers whose codes are the same, have an
// eigenvalue of one against the data's, which can be many orders of
// magnitude larger, and their error is as large as their residual. The
// coefficient matrix is the identity plus a positive semidefinite matrix,
// so the error of the effects is never larger than the residual, and the
// residual is therefore also required to be small next to the effects
// themselves: 1e-8 of their norm bounds the relative error of the effects
// by 1e-8. On the whitened scale every effect has the prior's unit
// variance, so the traits are on the same footing in that test.
//
// A direction can sit at rounding level above its tolerance, as one whose
// right hand side is zero does when it is coupled to the others: its
// tolerance is then a fraction of nothing, and the rounding it picks up
// from the larger directions cannot be removed. The tolerance is not
// floored at a fraction of the largest direction instead, because that
// would release a direction far smaller than the largest, such as a trait
// whose error variance is far larger than another's, long before it had
// converged. A restart only follows a recurrence that has converged, so
// when a restart fails to halve the true residual of every direction still
// above its tolerance, and that residual is at rounding level of the
// product it was computed from, the iteration stops, because more of it
// cannot help. A residual that has stopped falling while still above
// rounding level is not progress lost to rounding, and never stops it.
//
// Stopping there does not make the effects accurate when the variance
// ratio is so large that rounding moves the directions only the prior
// holds. The residual bounds their error, but loosely: rounding in the
// directions the records see dominates it and barely affects the effects.
// The restarts measure the error more directly. Each one solves for a
// correction to the effects, which is small when the residual is rounding
// in well determined directions and as large as the drift when the effects
// are still moving along poorly determined ones. A last correction above
// 1e-6 of the effects is reported in roundingBound rather than passed over.
//
// A restart that does not change the effects at all, while the residual is
// still above rounding level, has stalled: the preconditioner cannot
// resolve a direction the residual lies in, as when terms' covariance
// matrices are aligned to within the precision of its blocks. Further
// restarts would only repeat it, so the iteration stops and says so in
// stalledAt.
template<class ApplyCoef, class Precondition>
PcgStatus restartedPCG(arma::mat& At, const arma::mat& Bt,
                       ApplyCoef applyCoef, Precondition precondition,
                       double lamMax, arma::uword maxIter){
  const arma::uword nR = Bt.n_rows;
  const arma::uword nLoci = Bt.n_cols;
  arma::mat Rt = Bt;
  double rNorm0 = arma::norm(Bt, "fro");
  arma::vec rowTol = sqrt(sum(square(Bt),1));
  rowTol.elem(find(rowTol<arma::datum::eps*rNorm0)).fill(
    arma::datum::eps*rNorm0);
  rowTol *= 1.0e-8;
  const double errTol = 1.0e-8;
  At.zeros(nR,nLoci);
  auto rowsConverged = [&](const arma::mat& R){
    return arma::all(sqrt(sum(square(R),1))<=rowTol) &&
      (arma::norm(R, "fro")<=errTol*arma::norm(At, "fro"));
  };
  const arma::uword maxRestart = 20;
  arma::mat Zt(nR,nLoci), Apt(nR,nLoci), Pt;
  PcgStatus status{(rNorm0<=0.0), 0.0, 0.0};
  bool& converged = status.converged;
  arma::uword iter = 0;
  arma::vec prevNorm;
  double prevErr = arma::datum::inf;
  double lastMove = 0.0;
  for(arma::uword restart=0; (restart<=maxRestart)&&(!converged); ++restart){
    if(restart>0){
      // The true residual
      applyCoef(Apt, At);
      ++iter;
      Rt = Bt-Apt;
    }
    if(rowsConverged(Rt)){
      converged = true;
      break;
    }
    if(restart>0){
      arma::vec rowNorm = sqrt(sum(square(Rt),1));
      double rNorm = arma::norm(Rt, "fro");
      bool errOpen = rNorm>errTol*arma::norm(At, "fro");
      if(restart>1){
        arma::uvec open = find(rowNorm>rowTol);
        double scale = lamMax*arma::norm(At, "fro")+arma::norm(Bt, "fro");
        bool atRounding = rNorm<=64.0*arma::datum::eps*scale;
        if( !atRounding &&
            (lastMove<=arma::datum::eps*arma::norm(At, "fro")) ){
          status.stalledAt = rNorm/scale;
          break;
        }
        if(atRounding && arma::all(rowNorm(open)>0.5*prevNorm(open)) &&
           (!errOpen || (rNorm>0.5*prevErr))){
          converged = true;
          double move = lastMove/arma::norm(At, "fro");
          if(errOpen && (move>1.0e-6)){
            status.roundingBound = move;
          }
          break;
        }
      }
      prevNorm = rowNorm;
      prevErr = rNorm;
    }
    if(iter>=maxIter){
      break;
    }
    precondition(Zt, Rt);
    Pt = Zt;
    double rz = arma::accu(Rt%Zt);
    if(!(rz>0.0)){
      // The residual is at rounding level in every direction the
      // preconditioner sees
      converged = true;
      break;
    }
    bool claimed = false;
    arma::mat Aprev = At;
    while(iter<maxIter){
      applyCoef(Apt, Pt);
      ++iter;
      double pAp = arma::accu(Pt%Apt);
      if(!(pAp>0.0)){
        // The matrix is positive definite, so this is rounding once there is
        // nothing left to solve for
        claimed = true;
        break;
      }
      double step = rz/pAp;
      At += step*Pt;
      Rt -= step*Apt;
      if(rowsConverged(Rt)){
        claimed = true;
        break;
      }
      precondition(Zt, Rt);
      double rzNew = arma::accu(Rt%Zt);
      if(!(rzNew>0.0)){
        // Nothing the preconditioner can see is left in the residual, so
        // there is no direction left to search in
        claimed = true;
        break;
      }
      Pt = Zt+(rzNew/rz)*Pt;
      rz = rzNew;
    }
    lastMove = arma::norm(At-Aprev, "fro");
    if(!claimed){
      break;
    }
  }
  return status;
}

// Fits the multikernel model in fast mode. kinds lists the terms and
// domCovariate adds the summed dominance code to the fixed effects, as the
// direct form does for the models with dominance. VuIn, VeIn, estVarComp and
// subset are as for fastMKVarComp. nPrecond and maxRank set the deflation of
// the preconditioner, see mkDeflation; either of them zero leaves the block
// preconditioner alone.
FastMKFit fastMultiKernel(const arma::mat& Y, arma::uvec x,
                          const PackedGeno& g, const LociPos& pos,
                          const std::vector<TermKind>& kinds,
                          bool domCovariate,
                          const arma::field<arma::mat>& VuIn,
                          const arma::mat& VeIn, arma::uword maxIter,
                          bool estVarComp, arma::uvec subset,
                          arma::uword nPrecond, arma::uword maxRank,
                          int nThreads, const std::string& name){
  const arma::uword nT = Y.n_cols;
  const arma::uword nTerm = kinds.size();
  arma::uword ploidy = g.ploidy;
  arma::uword nInd = g.nInd;
  arma::uword nLoci = pos.base.size();
  if(Y.n_rows!=nInd){
    Rcpp::stop("%s: the responses must have one row per individual",
               name.c_str());
  }
  if(!Y.is_finite()){
    Rcpp::stop("%s: every trait must be recorded on every individual",
               name.c_str());
  }
  std::vector<FastTerm> terms;
  for(TermKind kind : kinds){
    terms.push_back(makeTerm(kind, ploidy));
  }
  arma::Mat<unsigned char> planes = packChipPlanes(g, pos, nThreads);

  // Fixed effects: the sum to zero design with an intercept, as in the
  // direct form, and for the models with dominance the summed dominance code
  arma::mat X = makeX(x);
  if(domCovariate){
    double scale = makeTerm(TermKind::dominance, ploidy).scale;
    X = join_rows(X, dominanceSum(planes, ploidy, scale, nThreads));
  }

  FastMKFit fit;
  fastMKVarComp(fit, Y, X, planes, ploidy, nLoci, terms, VuIn, VeIn,
                estVarComp, subset, nThreads, name);
  const Whitened w = whitenTerms(fit.Vu, fit.Ve, name);
  const std::vector<arma::mat>& F = w.F;
  const std::vector<arma::uword>& off = w.off;
  // With every variance zero there are no effects to solve for, nR is zero
  // and the model is fitted on its fixed effects alone, as the direct form
  // fits it
  const arma::uword nR = off[nTerm];

  // The fixed effects are absorbed with S = I-Q*Q', for an orthonormal basis
  // Q of their span. Columns of X that depend on the others, as the summed
  // dominance code can on the intercept, are dropped, as the direct form
  // drops them.
  FixedDesign fd = fixedDesign(X);
  const arma::mat& Q = fd.Q;

  fit.stats = markerStats(planes, ploidy, nLoci, terms, Q, nThreads);
  arma::mat offset = fit.stats.Mmean;
  for(arma::uword t=0; t<nTerm; ++t){
    offset.col(t) += terms[t].constant;
  }

  const BlockPrecond P = buildPrecond(fit.stats, w, nLoci, nThreads);
  bool anyNull = false;
  const std::vector<arma::mat> nullBasis =
    unseenDirections(fit.stats, w, nT, nInd, nLoci, nThreads, anyNull);
  const MKDeflation defl = mkDeflation(planes, ploidy, terms, w, P, nT, nInd,
                                       nLoci, nPrecond, maxRank, nThreads);
  fit.rank = defl.W.n_cols;

  arma::mat tabs(32*nT*nTerm, planes.n_cols);
  std::vector<arma::mat> Ut(nTerm), Wt(nTerm);
  for(arma::uword t=0; t<nTerm; ++t){
    Ut[t].set_size(nT,nLoci);
    Wt[t].set_size(nT,nLoci);
  }
  arma::mat Tt(nT,nInd);

  // One product with the coefficient matrix, out = (I+F'*W'*S*W*F)*in, with
  // the effects of every term stacked in the rows of in and out
  auto applyCoef = [&](arma::mat& out, const arma::mat& in){
    for(arma::uword t=0; t<nTerm; ++t){
      if(off[t+1]==off[t]){
        Ut[t].zeros();
      }else{
        Ut[t] = F[t]*in.rows(off[t], off[t+1]-1);
      }
    }
    multTerms(Tt, planes, ploidy, terms, offset, Ut, tabs, nThreads);
    Tt -= (Tt*Q)*Q.t(); // S applied to each trait
    multTermsT(Wt, planes, ploidy, terms, offset, Tt, nThreads);
    out = in;
    for(arma::uword t=0; t<nTerm; ++t){
      if(off[t+1]>off[t]){
        out.rows(off[t], off[t+1]-1) += F[t].t()*Wt[t];
      }
    }
    deflateUnseen(out, nullBasis, anyNull);
  };
  auto precondition = [&](arma::mat& Z, const arma::mat& R){
    if(defl.W.n_cols==0){
      applyPrecond(Z, R, P, nThreads);
    }else{
      applyDeflatedPrecond(Z, R, defl, nThreads);
    }
    deflateUnseen(Z, nullBasis, anyNull);
  };

  // Right hand side, F'*W'*S*Y* for the whitened responses Y* = Y*Le^-T,
  // from the responses with each level's offset taken off, with the fixed
  // effects absorbed before the result is whitened
  arma::mat Yc;
  const arma::mat levelOffset = levelOffsets(Y, x, Yc);
  const arma::uword nLevel = levelOffset.n_rows;
  Tt = w.Linv*(Yc-Q*(Q.t()*Yc)).t();
  multTermsT(Wt, planes, ploidy, terms, offset, Tt, nThreads);
  arma::mat Bt(nR,nLoci);
  for(arma::uword t=0; t<nTerm; ++t){
    if(off[t+1]>off[t]){
      Bt.rows(off[t], off[t+1]-1) = F[t].t()*Wt[t];
    }
  }
  deflateUnseen(Bt, nullBasis, anyNull);

  // The largest diagonal entry of the coefficient matrix, a lower bound on
  // its norm, for judging when a residual is at rounding level. Rounding in
  // a product with the coefficient matrix grows with the matrix times the
  // effects, not with the product, which is far smaller when the effects
  // cancel.
  double lamMax = 1.0;
  if(nR>0){
    lamMax = 1.0/square(P.ds).min();
  }
  arma::mat At;
  const PcgStatus status = restartedPCG(At, Bt, applyCoef, precondition,
                                        lamMax, maxIter);
  if(status.stalledAt>0.0){
    Rf_warning("%s: the iteration stalled without converging, with a "
               "residual of %.2g of the right hand sides", name.c_str(),
               status.stalledAt);
  }else if(!status.converged){
    Rf_warning("%s: reached maxIter without converging", name.c_str());
  }else if(status.roundingBound>0.0){
    Rf_warning("%s: rounding limits the accuracy of the solution; the "
               "marker effects still change by about %.2g of their size "
               "between restarts",
               name.c_str(), status.roundingBound);
  }

  deflateUnseen(At, nullBasis, anyNull);

  // Back to the marker effects of the original traits. On the whitened
  // scale the effects are F_k*A_k, held with one row per trait, and
  // U_k = U*_k*Le' undoes the whitening.
  fit.U.resize(nTerm);
  for(arma::uword t=0; t<nTerm; ++t){
    if(off[t+1]==off[t]){
      Ut[t].zeros();
    }else{
      Ut[t] = w.Le*(F[t]*At.rows(off[t], off[t+1]-1));
    }
    fit.U[t] = Ut[t].t();
  }

  // Fixed effects, recovered from the solution. With the same fixed effects
  // for every trait and every trait recorded, generalized least squares
  // reduces to least squares trait by trait. They are found against the
  // centered codes, while the direct form, and the genetic values, use the
  // uncentered ones, so the intercept takes off the means times the effects.
  multTerms(Tt, planes, ploidy, terms, offset, Ut, tabs, nThreads);
  // From the responses with each level's offset taken off, as for the right
  // hand side, and those offsets then added back. In the sum to zero design
  // of makeX the intercept is the mean of the levels and the first nLevel-1
  // contrasts are each level's departure from it, so the offsets add their
  // mean to the intercept and their departures from it to the contrasts.
  arma::mat Bq = Q.t()*(Yc-Tt.t());
  fit.beta = fixedEffectsFromBasis(fd, Bq, X.n_cols);
  if(nLevel>0){
    arma::rowvec offsetMean = arma::mean(levelOffset, 0);
    fit.beta.row(0) += offsetMean;
    for(arma::uword l=0; l+1<nLevel; ++l){
      fit.beta.row(1+l) += levelOffset.row(l)-offsetMean;
    }
  }
  for(arma::uword t=0; t<nTerm; ++t){
    fit.beta.row(0) -= fit.stats.Mmean.col(t).t()*fit.U[t];
  }
  return fit;
}

// The variance components in the form the R functions expect: for one trait,
// one row per term, and otherwise a list of covariance matrices
Rcpp::RObject vuOut(const FastMKFit& fit){
  arma::uword nTerm = fit.Vu.size();
  if(fit.Ve.n_rows==1){
    arma::mat Vu(nTerm,1);
    for(arma::uword t=0; t<nTerm; ++t){
      Vu(t,0) = fit.Vu[t](0,0);
    }
    return Rcpp::wrap(Vu);
  }
  Rcpp::List Vu(nTerm);
  for(arma::uword t=0; t<nTerm; ++t){
    Vu[t] = Rcpp::wrap(fit.Vu[t]);
  }
  return Vu;
}

// Markers at which every record has the same dosage, which the direct form
// gives no effects
arma::uvec fixedMarkers(const arma::mat& genoCount, double nInd){
  arma::uvec fixed(genoCount.n_cols);
  for(arma::uword j=0; j<genoCount.n_cols; ++j){
    fixed(j) = arma::any(genoCount.col(j)==nInd) ? 1 : 0;
  }
  return fixed;
}

} // namespace

// Called by RRBLUP_D with fastMode TRUE
// x is an indicator vector for the fixed effect levels, as in callRRBLUP
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP_D(const arma::mat& Y, arma::uvec x,
                            const Rcpp::S4& pop,
                            arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                            const arma::field<arma::mat>& Vu,
                            const arma::mat& Ve, arma::uword maxIter,
                            bool estVarComp, arma::uvec subset,
                            arma::uword nPrecond, arma::uword maxRank,
                            int nThreads){
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::additive, TermKind::dominance};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, true, Vu, Ve,
                                  maxIter, estVarComp, subset, nPrecond,
                                  maxRank, nThreads, "RRBLUP_D fastMode");
  // Average effects from the genotype frequencies, as in callRRBLUP_D
  double dn = double(g.nInd);
  arma::uvec fixed = fixedMarkers(fit.stats.genoCount, dn);
  arma::mat freqMat = fit.stats.genoCount/dn;
  arma::mat alpha, d;
  domAverageEffects(freqMat, fixed, fit.U[0], fit.U[1],
                    fit.beta.row(fit.beta.n_rows-1), g.ploidy, nThreads,
                    alpha, d);
  arma::rowvec Mmean = fit.stats.Mmean.col(0).t();
  return Rcpp::List::create(
    Rcpp::Named("alpha")=alpha,
    Rcpp::Named("beta")=negMeanProduct(Mmean,alpha),
    Rcpp::Named("a")=fit.U[0],
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=fit.beta.row(0),
    Rcpp::Named("Vu")=vuOut(fit),
    Rcpp::Named("Ve")=fit.Ve,
    Rcpp::Named("rank")=int(fit.rank)
  );
}

// Called by RRBLUP_GCA with fastMode TRUE
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP_GCA(const arma::mat& Y, arma::uvec x,
                              const Rcpp::S4& pop,
                              arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                              const arma::field<arma::mat>& Vu,
                              const arma::mat& Ve, arma::uword maxIter,
                              bool estVarComp, arma::uvec subset,
                              arma::uword nPrecond, arma::uword maxRank,
                              int nThreads){
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::maternal, TermKind::paternal};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, false, Vu, Ve,
                                  maxIter, estVarComp, subset, nPrecond,
                                  maxRank, nThreads, "RRBLUP_GCA fastMode");
  arma::rowvec Mmean1 = fit.stats.Mmean.col(0).t();
  arma::rowvec Mmean2 = fit.stats.Mmean.col(1).t();
  return Rcpp::List::create(Rcpp::Named("alpha1")=fit.U[0],
                            Rcpp::Named("alpha2")=fit.U[1],
                            Rcpp::Named("beta1")=negMeanProduct(Mmean1,fit.U[0]),
                            Rcpp::Named("beta2")=negMeanProduct(Mmean2,fit.U[1]),
                            Rcpp::Named("mu")=fit.beta.row(0),
                            Rcpp::Named("Vu")=vuOut(fit),
                            Rcpp::Named("Ve")=fit.Ve,
                            Rcpp::Named("rank")=int(fit.rank));
}

// Called by RRBLUP_SCA with fastMode TRUE
// Only works with diploids
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP_SCA(const arma::mat& Y, arma::uvec x,
                              const Rcpp::S4& pop,
                              arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                              const arma::field<arma::mat>& Vu,
                              const arma::mat& Ve, arma::uword maxIter,
                              bool estVarComp, arma::uvec subset,
                              arma::uword nPrecond, arma::uword maxRank,
                              int nThreads){
  PackedGeno g(pop);
  if(g.ploidy!=2){
    Rcpp::stop("RRBLUP_SCA fastMode: only diploids are supported");
  }
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::maternal, TermKind::paternal,
                                 TermKind::dominance};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, true, Vu, Ve,
                                  maxIter, estVarComp, subset, nPrecond,
                                  maxRank, nThreads, "RRBLUP_SCA fastMode");
  // Average effects of each parent, as in callRRBLUP_SCA
  double dn = double(g.nInd);
  arma::uvec fixed = fixedMarkers(fit.stats.genoCount, dn);
  arma::rowvec Mmean1 = fit.stats.Mmean.col(0).t();
  arma::rowvec Mmean2 = fit.stats.Mmean.col(1).t();
  arma::rowvec p12 = fit.stats.both11/dn;
  arma::mat alpha1, alpha2, d;
  scaAverageEffects(Mmean1, Mmean2, p12, fixed, fit.U[0], fit.U[1], fit.U[2],
                    fit.beta.row(fit.beta.n_rows-1), nThreads,
                    alpha1, alpha2, d);
  return Rcpp::List::create(
    Rcpp::Named("alpha1")=alpha1,
    Rcpp::Named("alpha2")=alpha2,
    Rcpp::Named("beta1")=negMeanProduct(Mmean1,alpha1),
    Rcpp::Named("beta2")=negMeanProduct(Mmean2,alpha2),
    Rcpp::Named("a1")=fit.U[0],
    Rcpp::Named("a2")=fit.U[1],
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=fit.beta.row(0),
    Rcpp::Named("Vu")=vuOut(fit),
    Rcpp::Named("Ve")=fit.Ve,
    Rcpp::Named("rank")=int(fit.rank)
  );
}
