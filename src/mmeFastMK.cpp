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

#include <algorithm>

#include "mme.h"

// Defined in mmeMultiKernel.cpp. Declared here without its default, which
// Rcpp reads from the definition.
Rcpp::List solveRRBLUPMK(const arma::mat& y, const arma::mat& X,
                         const arma::field<arma::mat>& Mlist,
                         int maxIter);

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
      for(arma::uword i=0; i<nInd; ++i){
        const double* ti = Tt.colptr(i);
        const unsigned char* row = col+i*ploidy;
        for(arma::uword t=0; t<nTerm; ++t){
          const FastTerm& term = terms[t];
          double* bt = bins.data()+(2*t+(i & 1u))*setLen;
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
// subset these are the codes themselves. Otherwise the cross product of the
// codes, M*M', is built a panel of markers at a time and replaced by a
// factor F with F*F' = M*M' over its positive eigenvalues, which has no more
// columns than there are records. The model depends on a term only through
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
  // eigen2 writes into outputs of the right size, as in remlFromGram
  arma::vec eigval(nSub);
  arma::mat U(nSub,nSub);
  eigen2(eigval, U, K);
  K.reset();
  arma::uvec keep = find(eigval>rankCutoff(eigval.max(), nSub));
  out = U.cols(keep);
  out.each_row() %= sqrt(eigval(keep)).t();
  return out;
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
};

// Fits the multikernel model in fast mode. kinds lists the terms and
// domCovariate adds the summed dominance code to the fixed effects, as the
// direct form does for the models with dominance. VuIn and VeIn are used
// when estVarComp is false, as nT by nT matrices with one for each term.
// Otherwise they are estimated by REML on the records in subset, which is
// one based.
FastMKFit fastMultiKernel(const arma::mat& Y, arma::uvec x,
                          const PackedGeno& g, const LociPos& pos,
                          const std::vector<TermKind>& kinds,
                          bool domCovariate,
                          const arma::field<arma::mat>& VuIn,
                          const arma::mat& VeIn, arma::uword maxIter,
                          bool estVarComp, arma::uvec subset, int nThreads,
                          const std::string& name){
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

  // Variance components, by REML on a subset of the records when they have
  // not been supplied. As in the fast mode of RRBLUP, the decomposition REML
  // needs costs the square of the number of records in memory and their cube
  // in time, so it is done on a subset, and the marker effects are then
  // solved for on every record.
  FastMKFit fit;
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
    Rcpp::List varComp = solveRRBLUPMK(ySub, XSub, Mlist, fastMKRemlIter);
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

  // Whitening by Ve
  arma::mat Le;
  if(!arma::chol(Le, fit.Ve, "lower")){
    Rcpp::stop("%s: Ve must be positive definite", name.c_str());
  }
  arma::mat Linv;
  if(!arma::inv(Linv, arma::trimatl(Le))){
    Rcpp::stop("%s: Ve must be positive definite", name.c_str());
  }
  // Each term's covariance on the whitened scale as F*F'. These are
  // eigenvalues of a covariance matrix among traits, so rankCutoff does not
  // apply, and a positive eigenvalue is kept however small it is next to
  // the others, for the reason given in callFastRRBLUP_MV. A tiny one only
  // adds a direction that barely moves.
  std::vector<arma::mat> F(nTerm);
  std::vector<arma::uword> off(nTerm+1, 0);
  for(arma::uword t=0; t<nTerm; ++t){
    arma::mat A = Linv*fit.Vu[t]*Linv.t();
    A = 0.5*(A+A.t());
    arma::vec lam;
    arma::mat Qv;
    if(!arma::eig_sym(lam, Qv, A)){
      Rcpp::stop("%s: decomposition of Vu failed", name.c_str());
    }
    arma::uvec keep = find(lam>0.0);
    F[t] = Qv.cols(keep);
    F[t].each_row() %= sqrt(lam(keep)).t();
    off[t+1] = off[t]+keep.n_elem;
  }
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

  // The preconditioner is block Jacobi, with one block for each marker over
  // every term and whitened trait. The terms read at one marker are often
  // correlated, as the additive and dominance codes are wherever the allele
  // frequency is not one half, so the blocks hold their cross products as
  // well as the diagonal. Each block is the identity plus a positive
  // semidefinite matrix, so it always has an inverse.
  arma::cube Binv(nR,nR,nLoci);
  if(nR>0){
    arma::uword nBlocks = countBlocks(nLoci);
    int precondThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(precondThreads)){
      precondThreads = static_cast<int>(nBlocks);
    }
    // F_k'*F_l for every pair of terms, the same for every marker
    std::vector<arma::mat> FF(nTerm*nTerm);
    for(arma::uword t=0; t<nTerm; ++t){
      for(arma::uword t2=0; t2<nTerm; ++t2){
        FF[t*nTerm+t2] = F[t].t()*F[t2];
      }
    }
    // Each block records its own failure, so no flag is shared by threads
    std::vector<int> failed(nBlocks, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(precondThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      arma::mat B(nR,nR), Bi;
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        B.eye();
        for(arma::uword t=0; t<nTerm; ++t){
          for(arma::uword t2=0; t2<nTerm; ++t2){
            if( (off[t+1]==off[t]) || (off[t2+1]==off[t2]) ){
              continue;
            }
            B.submat(off[t], off[t2], off[t+1]-1, off[t2+1]-1) +=
              fit.stats.Dkl(t*nTerm+t2,j)*FF[t*nTerm+t2];
          }
        }
        B = 0.5*(B+B.t());
        if(arma::inv_sympd(Bi, B)){
          Binv.slice(j) = Bi;
        }else{
          failed[block] = 1;
        }
      }
    }
    if(std::find(failed.begin(), failed.end(), 1)!=failed.end()){
      Rcpp::stop("%s: preconditioner could not be formed", name.c_str());
    }
  }

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
  };
  // Each marker's block is applied to its own column, so the work splits
  // over blocks of markers
  arma::uword precondBlocks = countBlocks(nLoci);
  int precondThreads = nThreads;
  if(precondBlocks < static_cast<arma::uword>(precondThreads)){
    precondThreads = static_cast<int>(precondBlocks);
  }
  auto precondition = [&](arma::mat& Z, const arma::mat& R){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(precondThreads)
#endif
    for(arma::uword block=0; block<precondBlocks; ++block){
      arma::uword lociStart = blockStart(nLoci, precondBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, precondBlocks, block+1);
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        const double* Bj = Binv.slice_memptr(j);
        const double* rj = R.colptr(j);
        double* zj = Z.colptr(j);
        for(arma::uword a=0; a<nR; ++a){
          double sum = 0.0;
          for(arma::uword c=0; c<nR; ++c){
            sum += Bj[a+c*nR]*rj[c];
          }
          zj[a] = sum;
        }
      }
    }
  };

  // Right hand side, F'*W'*S*Y* for the whitened responses Y* = Y*Le^-T
  Tt = Linv*Y.t();
  Tt -= (Tt*Q)*Q.t();
  multTermsT(Wt, planes, ploidy, terms, offset, Tt, nThreads);
  arma::mat Rt(nR,nLoci);
  for(arma::uword t=0; t<nTerm; ++t){
    if(off[t+1]>off[t]){
      Rt.rows(off[t], off[t+1]-1) = F[t].t()*Wt[t];
    }
  }

  // Preconditioned conjugate gradient on the stacked whitened effects. The
  // terms and traits share one recurrence, because the terms are coupled
  // through the records.
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
  // A direction can sit at rounding level above its tolerance, as one whose
  // right hand side is zero does when it is coupled to the others: its
  // tolerance is then a fraction of nothing, and the rounding it picks up
  // from the larger directions cannot be removed. The tolerance is not
  // floored at a fraction of the largest direction instead, because that
  // would release a direction far smaller than the largest, such as a trait
  // whose error variance is far larger than another's, long before it had
  // converged. A restart only follows a recurrence that has converged, so
  // when a restart fails to halve the true residual of every direction still
  // above its tolerance, those directions are at rounding level and the
  // solution is accepted as converged.
  const arma::mat Bt = Rt;
  double rNorm0 = arma::norm(Bt, "fro");
  arma::vec rowTol = sqrt(sum(square(Bt),1));
  rowTol.elem(find(rowTol<arma::datum::eps*rNorm0)).fill(
    arma::datum::eps*rNorm0);
  rowTol *= 1.0e-8;
  auto rowsConverged = [&](const arma::mat& R){
    return arma::all(sqrt(sum(square(R),1))<=rowTol);
  };
  const arma::uword maxRestart = 20;
  arma::mat At(nR,nLoci,arma::fill::zeros);
  arma::mat Zt(nR,nLoci), Apt(nR,nLoci), Pt;
  bool converged = (rNorm0<=0.0);
  arma::uword iter = 0;
  arma::vec prevNorm;
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
      if(restart>1){
        arma::uvec open = find(rowNorm>rowTol);
        if(arma::all(rowNorm(open)>0.5*prevNorm(open))){
          converged = true;
          break;
        }
      }
      prevNorm = rowNorm;
    }
    if(iter>=maxIter){
      break;
    }
    precondition(Zt, Rt);
    Pt = Zt;
    double rz = arma::accu(Rt%Zt);
    bool claimed = false;
    while(iter<maxIter){
      applyCoef(Apt, Pt);
      ++iter;
      double pAp = arma::accu(Pt%Apt);
      if(pAp<=0.0){
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
      Pt = Zt+(rzNew/rz)*Pt;
      rz = rzNew;
    }
    if(!claimed){
      break;
    }
  }
  if(!converged){
    Rf_warning("%s: reached maxIter without converging", name.c_str());
  }

  // Back to the marker effects of the original traits. On the whitened
  // scale the effects are F_k*A_k, held with one row per trait, and
  // U_k = U*_k*Le' undoes the whitening.
  fit.U.resize(nTerm);
  for(arma::uword t=0; t<nTerm; ++t){
    if(off[t+1]==off[t]){
      Ut[t].zeros();
    }else{
      Ut[t] = Le*(F[t]*At.rows(off[t], off[t+1]-1));
    }
    fit.U[t] = Ut[t].t();
  }

  // Fixed effects, recovered from the solution. With the same fixed effects
  // for every trait and every trait recorded, generalized least squares
  // reduces to least squares trait by trait. They are found against the
  // centered codes, while the direct form, and the genetic values, use the
  // uncentered ones, so the intercept takes off the means times the effects.
  multTerms(Tt, planes, ploidy, terms, offset, Ut, tabs, nThreads);
  arma::mat Bq = Q.t()*(Y-Tt.t());
  fit.beta = fixedEffectsFromBasis(fd, Bq, X.n_cols);
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
                            int nThreads){
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::additive, TermKind::dominance};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, true, Vu, Ve,
                                  maxIter, estVarComp, subset, nThreads,
                                  "RRBLUP_D fastMode");
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
    Rcpp::Named("Ve")=fit.Ve
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
                              int nThreads){
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::maternal, TermKind::paternal};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, false, Vu, Ve,
                                  maxIter, estVarComp, subset, nThreads,
                                  "RRBLUP_GCA fastMode");
  arma::rowvec Mmean1 = fit.stats.Mmean.col(0).t();
  arma::rowvec Mmean2 = fit.stats.Mmean.col(1).t();
  return Rcpp::List::create(Rcpp::Named("alpha1")=fit.U[0],
                            Rcpp::Named("alpha2")=fit.U[1],
                            Rcpp::Named("beta1")=negMeanProduct(Mmean1,fit.U[0]),
                            Rcpp::Named("beta2")=negMeanProduct(Mmean2,fit.U[1]),
                            Rcpp::Named("mu")=fit.beta.row(0),
                            Rcpp::Named("Vu")=vuOut(fit),
                            Rcpp::Named("Ve")=fit.Ve);
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
                              int nThreads){
  PackedGeno g(pop);
  if(g.ploidy!=2){
    Rcpp::stop("RRBLUP_SCA fastMode: only diploids are supported");
  }
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  std::vector<TermKind> kinds = {TermKind::maternal, TermKind::paternal,
                                 TermKind::dominance};
  FastMKFit fit = fastMultiKernel(Y, x, g, pos, kinds, true, Vu, Ve,
                                  maxIter, estVarComp, subset, nThreads,
                                  "RRBLUP_SCA fastMode");
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
    Rcpp::Named("Ve")=fit.Ve
  );
}
