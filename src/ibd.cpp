#include "alphasimr.h"
#include "ibd.h"
#include <algorithm>
#include <vector>

namespace {

// One two-column matrix of a packed record, as pointers to its columns
struct Block{
  const int* src;
  const int* site;
  int n;
};

// Read-only view of a packed record. See ibd.h for the layout. It holds
// plain pointers, so it can be read inside a parallel region once it has
// been built and checked.
class PackedView{
public:
  int nChr;
  int ploidy;
  std::vector<Block> blocks;

  PackedView() : nChr(0), ploidy(0) {}

  PackedView(const int* data, R_xlen_t len){
    if((len<2) || (data[0]>=0) || (data[1]<1)){
      Rcpp::stop("Invalid packed recombination record");
    }
    nChr = -data[0];
    ploidy = data[1];
    R_xlen_t nBlock = R_xlen_t(nChr)*ploidy;
    if(len<(2+nBlock)){
      Rcpp::stop("Invalid packed recombination record");
    }
    blocks.resize(nBlock);
    R_xlen_t pos = 2+nBlock;
    for(R_xlen_t i=0; i<nBlock; ++i){
      int n = data[2+i];
      if((n<1) || ((pos+2*R_xlen_t(n))>len)){
        Rcpp::stop("Invalid packed recombination record");
      }
      blocks[i].src = data+pos;
      blocks[i].site = data+pos+n;
      blocks[i].n = n;
      pos += 2*R_xlen_t(n);
    }
    if(pos!=len){
      Rcpp::stop("Invalid packed recombination record");
    }
  }

  const Block& get(int chr, int copy) const{
    return blocks[std::size_t(chr)*ploidy+copy];
  }
};

// Views a packed record held by R. Founders' records are not packed and are
// refused here, as are missing entries.
PackedView viewR(SEXP x){
  if((TYPEOF(x)!=INTSXP) || (Rf_xlength(x)<2)){
    Rcpp::stop("Invalid packed recombination record");
  }
  return PackedView(INTEGER(x), Rf_xlength(x));
}

// Segments of one chromosome copy being assembled
struct Segs{
  std::vector<int> src;
  std::vector<int> site;
};

// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns last element if value is greater than values of x
// Set left to the smallest value of the interval to search
int intervalSearchInt(const int* x, int n, int value){
  int left = 0;

  // Check if crossover is before beginning
  if(x[left]>value){
    return -1;
  }
  int end = n-1;

  // Check if crossover is at or past end
  if(x[end]<=value){
    return end;
  }

  // Perform search
  int right = end;
  while((right-left)>1){ // Interval can be decreased
    int middle = (left + right) / 2;
    if (x[middle] == value){
      left = middle;

      // Check if at the end of the vector
      if(left<end){
        // Check for identical values to the right
        while(x[left+1]==value){
          left += 1;
          if(left==end){
            break;
          }
        }
      }
      break;
    } else if (x[middle]>value){
      right = middle;
    }else{
      left = middle;
    }
  }
  return left;
}

// Appends rows from to to (inclusive) of a parental block
void appendRows(Segs& X, const Block& par, int from, int to){
  for(int r=from; r<=to; ++r){
    X.src.push_back(par.src[r]);
    X.site.push_back(par.site[r]);
  }
}

// Creates the first section of the IBD matrix
void firstSection(Segs& X, const Block& par, int recEnd){
  // Determine last site before new recombination starts
  // Assumes recEnd >= 2 and row 1 of par is always 1
  int matEnd = intervalSearchInt(par.site, par.n, recEnd-1);
  appendRows(X, par, 0, matEnd);
}

// Appends a section running from recStart up to, but not including, recEnd.
// lastSection passes the end of the parental block as matEnd instead.
void nextSection(Segs& X, const Block& par, int recStart, int matEnd){
  // Determine row at or just before recombination
  int matStart = intervalSearchInt(par.site, par.n, recStart);
  if((matStart<0) || X.src.empty()){
    Rcpp::stop("Invalid recombination history");
  }

  // Check if previous IBD section matches section after recombination
  if(par.src[matStart] == X.src.back()){
    // Move to the next IBD section in par. Anything left to copy is
    // appended, and nothing is when matStart has passed matEnd.
    appendRows(X, par, matStart+1, matEnd);
    return;
  }

  // Determine if recombination matches start of IBD section (probably rare)
  // Cases where this happens without an IBD change are handled above
  if(par.site[matStart] == recStart){
    // Relevant rows of par can be appended directly
    appendRows(X, par, matStart, matEnd);
    return;
  }

  // Recombination occurred within an IBD block, creating a new segment
  X.src.push_back(par.src[matStart]);
  X.site.push_back(recStart);
  appendRows(X, par, matStart+1, matEnd);
}

void middleSection(Segs& X, const Block& par, int recStart, int recEnd){
  // Determine last site before new recombination starts
  // May equal the row at the start
  int matEnd = intervalSearchInt(par.site, par.n, recEnd-1);
  nextSection(X, par, recStart, matEnd);
}

void lastSection(Segs& X, const Block& par, int recStart){
  nextSection(X, par, recStart, par.n-1);
}

// Packs the segments of every chromosome copy, ordered by chromosome and
// then by copy, into the layout described in ibd.h
std::vector<int> packSegs(const std::vector<Segs>& segs, int nChr, int ploidy){
  std::size_t len = 2+segs.size();
  for(const Segs& s : segs){
    len += 2*s.src.size();
  }
  std::vector<int> out;
  out.reserve(len);
  out.push_back(-nChr);
  out.push_back(ploidy);
  for(const Segs& s : segs){
    out.push_back(int(s.src.size()));
  }
  for(const Segs& s : segs){
    out.insert(out.end(), s.src.begin(), s.src.end());
    out.insert(out.end(), s.site.begin(), s.site.end());
  }
  return out;
}

// IBD of a founder, whose record lists one founder haplotype per copy. Each
// copy descends from its own haplotype along the whole of every chromosome.
std::vector<int> founderIbd(const Rcpp::IntegerVector& founder, int nChr){
  int ploidy = founder.size();
  std::vector<Segs> segs(std::size_t(nChr)*ploidy);
  for(int i=0; i<nChr; ++i){
    for(int k=0; k<ploidy; ++k){
      Segs& s = segs[std::size_t(i)*ploidy+k];
      s.src.push_back(founder[k]);
      s.site.push_back(1);
    }
  }
  return packSegs(segs, nChr, ploidy);
}

// IBD of a non-founder from its recombination history and its parents' IBD
std::vector<int> nonFounderIbd(const PackedView& rec,
                               const PackedView& mother,
                               const PackedView& father){
  int nChr = rec.nChr;
  int ploidy = rec.ploidy;
  int motherPloidy = mother.ploidy;
  int fatherPloidy = father.ploidy;
  if((mother.nChr!=nChr) || (father.nChr!=nChr)){
    Rcpp::stop("Parents have a different number of chromosomes");
  }

  // Test ploidy levels of parents
  if(ploidy == (motherPloidy+fatherPloidy)/2){
    // Suspect regular cross or DH
    // Proceed as normal
  }else if(ploidy == (motherPloidy+fatherPloidy)/4){
    // Suspect reduceGenome was used
    // Proceed as normal, all gametes will be taken from mother
  }else if(ploidy == (motherPloidy+fatherPloidy)){
    // Suspect doubleGenome or mergeGenome was used
    // Double motherPloidy to pull all gametes
    motherPloidy *= 2;
  }else{
    // No idea what happened
    Rcpp::stop("Unexpected parental ploidy levels");
  }

  std::vector<Segs> segs(std::size_t(nChr)*ploidy);
  for(int i=0; i<nChr; ++i){
    for(int j=0; j<ploidy; ++j){
      const PackedView& par = (j<motherPloidy/2) ? mother : father;
      const Block& h = rec.get(i, j);
      for(int k=0; k<h.n; ++k){
        if((h.src[k]<1) || (h.src[k]>par.ploidy)){
          Rcpp::stop("Recombination history names a missing parental chromosome");
        }
      }
      Segs& X = segs[std::size_t(i)*ploidy+j];
      int nCO = h.n-1;
      if(nCO==0){
        // Direct copy of chromosome IBD
        const Block& p = par.get(i, h.src[0]-1);
        appendRows(X, p, 0, p.n-1);
      }else{
        // Transfer first section
        firstSection(X, par.get(i, h.src[0]-1), h.site[1]);

        // Transfer middle section(s)
        for(int k=1; k<nCO; ++k){
          middleSection(X, par.get(i, h.src[k]-1), h.site[k], h.site[k+1]);
        }

        // Transfer last section
        lastSection(X, par.get(i, h.src[nCO]-1), h.site[nCO]);
      }
    }
  }
  return packSegs(segs, nChr, ploidy);
}

} // namespace

Rcpp::IntegerVector packRecord(
    const arma::field<arma::field<arma::Mat<int> > >& x){
  int nChr = x.n_elem;
  int ploidy = x(0).n_elem;
  R_xlen_t len = 2+R_xlen_t(nChr)*ploidy;
  for(int i=0; i<nChr; ++i){
    for(int j=0; j<ploidy; ++j){
      len += x(i)(j).n_elem;
    }
  }
  Rcpp::IntegerVector out(Rcpp::no_init(len));
  int* o = out.begin();
  *o++ = -nChr;
  *o++ = ploidy;
  for(int i=0; i<nChr; ++i){
    for(int j=0; j<ploidy; ++j){
      *o++ = x(i)(j).n_rows;
    }
  }
  for(int i=0; i<nChr; ++i){
    for(int j=0; j<ploidy; ++j){
      // Armadillo is column major, so the source column comes first
      o = std::copy(x(i)(j).begin(), x(i)(j).end(), o);
    }
  }
  return out;
}

// Packs recombination histories or IBD haplotypes given as nested lists of
// two-column matrices (individual, chromosome, copy). Entries that are not
// lists, such as founders' records, are returned unchanged.
// [[Rcpp::export]]
Rcpp::List packHist(const Rcpp::List& x){
  R_xlen_t n = x.size();
  Rcpp::List out(n);
  for(R_xlen_t i=0; i<n; ++i){
    SEXP e = x[i];
    if(TYPEOF(e)!=VECSXP){
      out[i] = e;
      continue;
    }
    Rcpp::List chrList(e);
    if(chrList.size()==0){
      Rcpp::stop("A recombination record has no chromosomes");
    }
    arma::field<arma::field<arma::Mat<int> > > rec(chrList.size());
    for(R_xlen_t j=0; j<chrList.size(); ++j){
      Rcpp::List copyList(chrList[j]);
      rec(j).set_size(copyList.size());
      for(R_xlen_t k=0; k<copyList.size(); ++k){
        rec(j)(k) = Rcpp::as<arma::Mat<int> >(copyList[k]);
        if((rec(j)(k).n_cols!=2) || (rec(j)(k).n_rows<1)){
          Rcpp::stop("A recombination record must be a two column matrix");
        }
      }
      if(rec(j).n_elem!=rec(0).n_elem){
        Rcpp::stop("A recombination record has inconsistent ploidy");
      }
    }
    if(rec(0).n_elem==0){
      Rcpp::stop("A recombination record has no chromosome copies");
    }
    out[i] = packRecord(rec);
  }
  return out;
}

// Reverses packHist, returning nested lists of two-column integer matrices.
// Entries that are not packed records, such as founders' records and
// missing entries, are returned unchanged.
// [[Rcpp::export]]
Rcpp::List unpackHist(const Rcpp::List& x){
  R_xlen_t n = x.size();
  Rcpp::List out(n);
  for(R_xlen_t i=0; i<n; ++i){
    SEXP e = x[i];
    if((TYPEOF(e)!=INTSXP) || (Rf_xlength(e)<1) || (INTEGER(e)[0]>=0)){
      out[i] = e;
      continue;
    }
    PackedView v = viewR(e);
    Rcpp::List chrList(v.nChr);
    for(int j=0; j<v.nChr; ++j){
      Rcpp::List copyList(v.ploidy);
      for(int k=0; k<v.ploidy; ++k){
        const Block& b = v.get(j, k);
        Rcpp::IntegerMatrix m(b.n, 2);
        std::copy(b.src, b.src+b.n, m.begin());
        std::copy(b.site, b.site+b.n, m.begin()+b.n);
        copyList[k] = m;
      }
      // The nested lists were once returned straight from an Armadillo
      // field, which gives each list a one column dim attribute. It is kept
      // so the expanded form is exactly what was returned before.
      copyList.attr("dim") = Rcpp::Dimension(v.ploidy, 1);
      chrList[j] = copyList;
    }
    chrList.attr("dim") = Rcpp::Dimension(v.nChr, 1);
    out[i] = chrList;
  }
  return out;
}

// Calculates IBD haplotypes for a set of individuals in one call.
//
// uid lists the individuals to compute, sorted by iid. An individual's
// parents always have smaller iids, so working in this order means every
// parent is ready before its progeny, either because it was computed
// earlier in this call or because it is already in hap, the cache of IBD
// haplotypes indexed by iid. recHist, mother, father and isFounder each
// hold one entry per element of uid. Founders' records list their founder
// haplotypes, and their parents are ignored.
//
// Returns the packed IBD haplotypes of uid, in the same order.
// [[Rcpp::export]]
Rcpp::List computeIbd(const Rcpp::IntegerVector& uid,
                      const Rcpp::IntegerVector& mother,
                      const Rcpp::IntegerVector& father,
                      const Rcpp::LogicalVector& isFounder,
                      const Rcpp::List& recHist,
                      const Rcpp::List& hap,
                      int nChr){
  R_xlen_t n = uid.size();
  if((mother.size()!=n) || (father.size()!=n) ||
     (isFounder.size()!=n) || (recHist.size()!=n)){
    Rcpp::stop("uid, mother, father, isFounder and recHist must match");
  }
  for(R_xlen_t i=1; i<n; ++i){
    if(uid[i]<=uid[i-1]){
      Rcpp::stop("uid must be sorted and unique");
    }
  }
  std::vector<std::vector<int> > result(n);

  // Finds a parent's IBD, from this call or else from the cache
  auto parentView = [&](int iid, R_xlen_t child) -> PackedView {
    const int* begin = uid.begin();
    const int* found = std::lower_bound(begin, begin+child, iid);
    if((found!=(begin+child)) && (*found==iid)){
      const std::vector<int>& r = result[found-begin];
      return PackedView(r.data(), R_xlen_t(r.size()));
    }
    if((iid<1) || (iid>hap.size())){
      Rcpp::stop("Missing IBD haplotypes for a parent");
    }
    return viewR(hap[iid-1]);
  };

  for(R_xlen_t i=0; i<n; ++i){
    if(isFounder[i]){
      result[i] = founderIbd(Rcpp::as<Rcpp::IntegerVector>(recHist[i]), nChr);
    }else{
      PackedView rec = viewR(recHist[i]);
      PackedView m = parentView(mother[i], i);
      PackedView f = parentView(father[i], i);
      result[i] = nonFounderIbd(rec, m, f);
    }
  }

  Rcpp::List out(n);
  for(R_xlen_t i=0; i<n; ++i){
    out[i] = Rcpp::IntegerVector(result[i].begin(), result[i].end());
  }
  return out;
}

// Expands packed IBD haplotypes into a matrix of founder haplotypes, with
// one row per haplotype and one column per locus
// [[Rcpp::export]]
arma::Mat<int> createIbdMat(const Rcpp::List& ibd,
                            arma::uvec chr,
                            arma::uvec nLoci,
                            arma::uword ploidy,
                            arma::uword nThreads){
  // R to C++
  if((chr.n_elem==0) || (chr.min()<1)){
    Rcpp::stop("Invalid chromosomes");
  }
  chr -= 1;
  arma::uword nChr = chr.n_elem;
  arma::uword nInd = ibd.size();
  arma::uword totLoci = accu(nLoci(chr));

  // The records are checked and viewed before the parallel region, which
  // then reads them only through plain pointers
  std::vector<PackedView> views(nInd);
  for(arma::uword i=0; i<nInd; ++i){
    views[i] = viewR(ibd[i]);
    if((arma::uword(views[i].ploidy)!=ploidy) ||
       (arma::uword(views[i].nChr)<=chr.max())){
      Rcpp::stop("IBD haplotypes do not match the population");
    }
  }

  arma::Mat<int> output(totLoci,nInd*ploidy);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<nInd; ++i){
    for(arma::uword j=0; j<ploidy; ++j){
      arma::uword stop,start=0;
      for(arma::uword k=0; k<nChr; ++k){
        const Block& b = views[i].get(chr(k), j);
        arma::uword nSeg = b.n;
        if(nSeg>1){
          // First segments
          for(arma::uword l=0; l<(nSeg-1); ++l){
            stop = start + b.site[l+1] - b.site[l] - 1;
            output.col(i*ploidy+j).rows(start,stop).fill(b.src[l]);
            start = stop + 1;
          }
        }
        // Last segment
        stop = start + nLoci(chr(k)) - b.site[nSeg-1];
        output.col(i*ploidy+j).rows(start,stop).fill(b.src[nSeg-1]);
        start = stop+1;
      }
    }
  }

  return output.t();
}
