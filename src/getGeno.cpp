#include "alphasimr.h"
#include <algorithm>

/*
 * Genotype data is stored in a list of raw arrays, one per chromosome.
 * Each array has dimensions nLoci/8 by ploidy by nInd
 * Output returned with dimensions nInd by nLoci
 *
 * The arrays are read in place through PackedGeno, which never copies them
 * and cannot write to them. See packedGv.h for why an arma::field is not
 * used.
 */

namespace {

// The packed array stores the loci of one haplotype next to each other,
// which is what meiosis wants, while the output matrix is column major and
// so stores the individuals of one locus next to each other. Whichever of
// the two is walked along its fast axis, the other is stepped through with a
// stride, and for a large population that stride is a cache line per locus.
// The copy is therefore done in tiles small enough to hold both ends in
// cache: kIndBlock individuals by kLociBlock loci.
const arma::uword kIndBlock = 64;
const arma::uword kLociBlock = 256;

// Individuals per tile. Falls below kIndBlock when there are too few
// individuals to give every thread a tile of the full size.
inline arma::uword indBlockSize(arma::uword nInd, int nThreads){
  arma::uword block = kIndBlock;
  if(nThreads>1){
    arma::uword nT = (arma::uword) nThreads;
    arma::uword perThread = (nInd+nT-1)/nT;
    if(perThread<block){
      block = (perThread>0) ? perThread : 1;
    }
  }
  return block;
}

// Writes the allele dosage carried on haplotypes [pStart,pEnd) of one tile
// into output(ind, outCol0+j). genoPtr is the start of a chromosome's
// array, which has nBins bytes per haplotype.
//
// Raw pointers are used throughout because Armadillo's element accessors are
// bounds checked in this package, and the check would be paid once per locus
// per haplotype. The byte holding a locus is re-read for each haplotype
// rather than cached across loci; within a tile it is already in L1, and
// re-reading is measurably faster than carrying the branch that a cache
// needs.
inline void gatherDosage(const unsigned char* genoPtr,
                         arma::uword nBins, arma::uword ploidy,
                         const arma::uvec& chrLociLoc,
                         arma::uword pStart, arma::uword pEnd,
                         arma::uword indStart, arma::uword indEnd,
                         arma::uword lociStart, arma::uword lociEnd,
                         arma::uword outCol0,
                         arma::Mat<unsigned char>& output){
  const arma::uword nOutRow = output.n_rows;
  const arma::uword* lociPtr = chrLociLoc.memptr();
  unsigned char* outPtr = output.memptr();
  for(arma::uword ind=indStart; ind<indEnd; ++ind){
    const unsigned char* hap = genoPtr + nBins*ploidy*ind;
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const arma::uword byte = lociPtr[j]/8;
      const arma::uword bit = lociPtr[j]%8;
      unsigned int dose = 0;
      for(arma::uword p=pStart; p<pEnd; ++p){
        dose += (unsigned int) ((hap[byte+nBins*p]>>bit) & 1u);
      }
      outPtr[ind+nOutRow*(outCol0+j)] = (unsigned char) dose;
    }
  }
}

// As gatherDosage, but keeps the haplotypes apart, writing the allele of
// haplotype p into output(ind*ploidy+p, outCol0+j).
inline void gatherHaplo(const unsigned char* genoPtr,
                        arma::uword nBins, arma::uword ploidy,
                        const arma::uvec& chrLociLoc,
                        arma::uword indStart, arma::uword indEnd,
                        arma::uword lociStart, arma::uword lociEnd,
                        arma::uword outCol0,
                        arma::Mat<unsigned char>& output){
  const arma::uword nOutRow = output.n_rows;
  const arma::uword* lociPtr = chrLociLoc.memptr();
  unsigned char* outPtr = output.memptr();
  for(arma::uword ind=indStart; ind<indEnd; ++ind){
    const unsigned char* hap = genoPtr + nBins*ploidy*ind;
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const arma::uword byte = lociPtr[j]/8;
      const arma::uword bit = lociPtr[j]%8;
      unsigned char* outCol = outPtr + nOutRow*(outCol0+j) + ind*ploidy;
      for(arma::uword p=0; p<ploidy; ++p){
        outCol[p] = (unsigned char) ((hap[byte+nBins*p]>>bit) & 1u);
      }
    }
  }
}

// Checks that a set of loci fits the genotypes it is about to be read
// from, given the number of bytes per haplotype of each chromosome. The
// gather functions above use raw pointers and do no bounds checking of
// their own, so a locus outside its chromosome would read or write memory
// that does not belong to the genotypes. The LociMap validity method cannot
// catch this, because it never sees the genotypes. The check runs on the
// calling thread, before any parallel region, because an exception must not
// leave an OpenMP region. lociLoc is still one based here, so a zero or
// negative locus from R has become 0 or a very large value and is caught by
// the same comparison.
void checkLoci(const std::vector<arma::uword>& nBins,
               const arma::Col<int>& lociPerChr,
               const arma::uvec& lociLoc){
  if(lociPerChr.n_elem!=nBins.size()){
    Rcpp::stop("Length of lociPerChr does not match the number of chromosomes");
  }
  arma::uword start = 0;
  for(arma::uword i=0; i<lociPerChr.n_elem; ++i){
    if(lociPerChr(i)<0){
      Rcpp::stop("lociPerChr contains a negative value");
    }
    arma::uword nLociChr = (arma::uword) lociPerChr(i);
    if(nLociChr>lociLoc.n_elem-start){
      Rcpp::stop("sum(lociPerChr) is greater than length(lociLoc)");
    }
    arma::uword maxLoc = nBins[i]*8;
    for(arma::uword j=start; j<start+nLociChr; ++j){
      if((lociLoc(j)<1) || (lociLoc(j)>maxLoc)){
        Rcpp::stop("A locus in lociLoc is outside its chromosome");
      }
    }
    start += nLociChr;
  }
  if(start!=lociLoc.n_elem){
    Rcpp::stop("sum(lociPerChr) is less than length(lociLoc)");
  }
}

// Number of copies of the counted allele at each of a set of loci, summed
// over every haplotype of every individual. Nothing the size of the
// genotype matrix is formed. The loci are split over a fixed number of
// blocks and each block counts only its own loci, so the counts do not
// depend on the number of threads. They are whole numbers in any case.
arma::uvec countAlleles(const PackedGeno& g,
                        const arma::Col<int>& lociPerChr,
                        const arma::uvec& lociLoc, int nThreads){
  checkLoci(g.nBins, lociPerChr, lociLoc);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  arma::uword nLoci = pos.base.size();
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  arma::uvec count(nLoci, arma::fill::zeros);
  arma::uword nBlocks = countBlocks(nLoci);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  if(nThreads < 1){
    nThreads = 1;
  }
  arma::uword* countPtr = count.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    // Individuals are the outer loop, so the loci of one individual, which
    // sit close together in memory, are read one after another
    for(arma::uword i=0; i<nInd; ++i){
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        countPtr[j] += dosageAt(pos.base[j]+pos.indStride(j)*i,
                                pos.hapStride(j), ploidy, pos.bit(j));
      }
    }
  }
  return count;
}

// Allele frequencies from allele counts, worked out exactly as the mean of
// a matrix of dosages divided by the ploidy level would be. The counts are
// whole numbers well below 2^53, so a sum of dosages in any order is the
// same count, and dividing by nInd and then by the ploidy level repeats
// Armadillo's mean, which is a sum divided by the number of rows.
arma::rowvec freqFromCounts(const arma::uvec& count, arma::uword nInd,
                            arma::uword ploidy){
  if(nInd==0){
    // Whatever Armadillo makes of the mean of no rows, so that an empty
    // population behaves as it did when a dosage matrix was averaged
    arma::rowvec output = arma::mean(arma::mat(0,count.n_elem),0)/ploidy;
    return output;
  }
  arma::rowvec output(count.n_elem);
  for(arma::uword j=0; j<count.n_elem; ++j){
    output(j) = (double(count(j))/double(nInd))/double(ploidy);
  }
  return output;
}

} // namespace

// Allele dosage carried on haplotypes [pStart,pEnd) at a set of loci, as an
// nInd by nLoci matrix. pStart=0 and pEnd=ploidy give the genotype dosage.
// lociLoc is one based, as it comes from R.
arma::Mat<unsigned char> genoDosage(const PackedGeno& g,
                                    const arma::Col<int>& lociPerChr,
                                    arma::uvec lociLoc,
                                    arma::uword pStart, arma::uword pEnd,
                                    int nThreads){
  checkLoci(g.nBins, lociPerChr, lociLoc);
  // R to C++ index correction
  lociLoc -= 1;

  arma::uword nInd = g.nInd;
  arma::uword nChr = g.chrPtr.size();
  arma::uword ploidy = g.ploidy;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::uword indBlock = indBlockSize(nInd,nThreads);
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr),arma::fill::zeros);
  // One parallel region covers every chromosome, so threads start once
  // per call rather than once per chromosome. Every thread walks the
  // chromosome loop and works out the same loci offsets, and the work
  // sharing construct below splits the individuals between them.
#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
#endif
  {
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      const unsigned char* chrPtr = g.chrPtr[i];
      arma::uword nBins = g.nBins[i];
      arma::uword nLociChr = chrLociLoc.n_elem;
      arma::uword nIndBlock = (nInd+indBlock-1)/indBlock;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(arma::uword ib=0; ib<nIndBlock; ++ib){
        arma::uword indStart = ib*indBlock;
        arma::uword indEnd = std::min(indStart+indBlock,nInd);
        for(arma::uword j0=0; j0<nLociChr; j0+=kLociBlock){
          arma::uword j1 = std::min(j0+kLociBlock,nLociChr);
          gatherDosage(chrPtr,nBins,ploidy,chrLociLoc,pStart,pEnd,
                       indStart,indEnd,j0,j1,
                       (arma::uword) loc1,output);
        }
      }
    }
  }
  }
  return output;
}

// Haplotypes at a set of loci, as an nInd*ploidy by nLoci matrix.
// lociLoc is one based, as it comes from R.
arma::Mat<unsigned char> genoHaplo(const PackedGeno& g,
                                   const arma::Col<int>& lociPerChr,
                                   arma::uvec lociLoc, int nThreads){
  checkLoci(g.nBins, lociPerChr, lociLoc);
  // R to C++ index correction
  lociLoc -= 1;

  arma::uword nInd = g.nInd;
  arma::uword nChr = g.chrPtr.size();
  arma::uword ploidy = g.ploidy;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::uword indBlock = indBlockSize(nInd,nThreads);
  arma::Mat<unsigned char> output(nInd*ploidy,arma::sum(lociPerChr));
  // One parallel region covers every chromosome, so threads start once
  // per call rather than once per chromosome. Every thread walks the
  // chromosome loop and works out the same loci offsets, and the work
  // sharing construct below splits the individuals between them.
#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
#endif
  {
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      const unsigned char* chrPtr = g.chrPtr[i];
      arma::uword nBins = g.nBins[i];
      arma::uword nLociChr = chrLociLoc.n_elem;
      arma::uword nIndBlock = (nInd+indBlock-1)/indBlock;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(arma::uword ib=0; ib<nIndBlock; ++ib){
        arma::uword indStart = ib*indBlock;
        arma::uword indEnd = std::min(indStart+indBlock,nInd);
        for(arma::uword j0=0; j0<nLociChr; j0+=kLociBlock){
          arma::uword j1 = std::min(j0+kLociBlock,nLociChr);
          gatherHaplo(chrPtr,nBins,ploidy,chrLociLoc,
                      indStart,indEnd,j0,j1,
                      (arma::uword) loc1,output);
        }
      }
    }
  }
  }
  return output;
}

// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const Rcpp::List& geno,
                                 const arma::Col<int>& lociPerChr,
                                 arma::uvec lociLoc, int nThreads){
  PackedGeno g(geno);
  return genoDosage(g, lociPerChr, lociLoc, 0, g.ploidy, nThreads);
}

// [[Rcpp::export]]
arma::Mat<unsigned char> getMaternalGeno(const Rcpp::List& geno,
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads){
  PackedGeno g(geno);
  return genoDosage(g, lociPerChr, lociLoc, 0, g.ploidy/2, nThreads);
}

// [[Rcpp::export]]
arma::Mat<unsigned char> getPaternalGeno(const Rcpp::List& geno,
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads){
  PackedGeno g(geno);
  return genoDosage(g, lociPerChr, lociLoc, g.ploidy/2, g.ploidy, nThreads);
}

// Returns haplotype data in a matrix of nInd*ploidy by nLoci
// [[Rcpp::export]]
arma::Mat<unsigned char> getHaplo(const Rcpp::List& geno,
                                  const arma::Col<int>& lociPerChr,
                                  arma::uvec lociLoc, int nThreads){
  PackedGeno g(geno);
  return genoHaplo(g, lociPerChr, lociLoc, nThreads);
}

// Returns haplotype data in a matrix of nInd by nLoci for a single
// chromosome group. i.e. just female or male chromosomes for diploids
// [[Rcpp::export]]
arma::Mat<unsigned char> getOneHaplo(const Rcpp::List& geno,
                                     const arma::Col<int>& lociPerChr,
                                     arma::uvec lociLoc, int haplo, int nThreads){
  PackedGeno g(geno);
  checkLoci(g.nBins, lociPerChr, lociLoc);
  // haplo picks the haplotype that gatherDosage reads through a raw
  // pointer, so it is checked here and not left to the R wrappers alone.
  // NA_integer_ is the most negative int and fails the first comparison.
  if((haplo<1) || (static_cast<arma::uword>(haplo)>g.ploidy)){
    Rcpp::stop("haplo must be between 1 and the ploidy level");
  }
  arma::uword hapIdx = (arma::uword) (haplo-1); // R to C++
  return genoDosage(g, lociPerChr, lociLoc, hapIdx, hapIdx+1, nThreads);
}

// Manually sets haplotype data
//
// geno is taken as an R list rather than as an Armadillo field, because
// RcppArmadillo builds the field's cubes on R's own memory instead of
// copying it. Writing into them would change the population passed in, and
// every R object sharing its genotypes, as well as the one returned. Each
// chromosome with loci to set is therefore duplicated before it is written
// to, and the others are shared with the input, which R copies on write.
// [[Rcpp::export]]
Rcpp::List setHaplo(const Rcpp::List& geno,
                    const arma::Mat<unsigned char>& haplo,
                    const arma::Col<int>& lociPerChr,
                    arma::uvec lociLoc, int nThreads){
  arma::uword nChr = geno.size();
  if(nChr==0){
    Rcpp::stop("geno must contain at least one chromosome");
  }
  // Every chromosome is written through raw pointers below, so each one is
  // checked to be a raw array with the same ploidy and individuals.
  arma::uvec nBins(nChr);
  arma::uword ploidy = 0, nInd = 0;
  for(arma::uword i=0; i<nChr; ++i){
    SEXP chr = geno[i];
    SEXP dims = Rf_getAttrib(chr, R_DimSymbol);
    if((TYPEOF(chr)!=RAWSXP) || (Rf_length(dims)!=3)){
      Rcpp::stop("Each chromosome in geno must be a three dimensional raw array");
    }
    const int* d = INTEGER(dims);
    if(i==0){
      ploidy = d[1];
      nInd = d[2];
    }else if(((arma::uword) d[1]!=ploidy) || ((arma::uword) d[2]!=nInd)){
      Rcpp::stop("Every chromosome in geno must have the same ploidy and individuals");
    }
    nBins(i) = d[0];
  }
  checkLoci(std::vector<arma::uword>(nBins.begin(), nBins.end()),
            lociPerChr, lociLoc);
  // R to C++ index correction
  lociLoc -= 1;

  // haplo is read with bounds checked accessors inside the parallel region
  // below, where the exception they throw cannot be caught safely, so its
  // shape is checked here first.
  if((haplo.n_rows!=nInd*ploidy) || (haplo.n_cols!=lociLoc.n_elem)){
    Rcpp::stop("haplo must have nInd*ploidy rows and one column per locus");
  }
  // The copies are made on this thread, because R's allocator must not be
  // called from inside a parallel region.
  Rcpp::List output(nChr);
  std::vector<unsigned char*> outPtr(nChr, nullptr);
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      Rcpp::RawVector chr = Rcpp::clone(Rcpp::RawVector(geno[i]));
      outPtr[i] = chr.begin();
      output[i] = chr;
    }else{
      output[i] = geno[i];
    }
  }
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  // One parallel region covers every chromosome, so threads start once
  // per call rather than once per chromosome. Every thread walks the
  // chromosome loop and works out the same loci offsets, and the work
  // sharing construct below splits the individuals between them.
#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
#endif
  {
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
      const arma::uword chrBins = nBins(i);
      unsigned char* chrPtr = outPtr[i];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=0; p<ploidy; ++p){
          // The cube is stored as bins by ploidy by individuals
          unsigned char* hap = chrPtr + chrBins*(p + ploidy*ind);
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(hap[currentByte]);
          workBits[chrLociLoc(0)%8] = haplo(ind*ploidy+p,loc1);
          hap[currentByte] = toByte(workBits);
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(hap[currentByte]);
            }
            workBits[chrLociLoc(j)%8] = haplo(ind*ploidy+p,j+loc1);
            hap[currentByte] = toByte(workBits);
          }
        }
      }
    }
  }
  }
  return output;
}


// [[Rcpp::export]]
void writeGeno(const Rcpp::List& geno,
               const arma::Col<int>& lociPerChr,
               arma::uvec lociLoc,
               Rcpp::String filePath, int nThreads){
  arma::Mat<unsigned char> output;
  output = getGeno(geno,lociPerChr,lociLoc,nThreads);
  std::ofstream outFile;
  outFile.open(filePath, std::ios_base::app);
  output.save(outFile,arma::raw_ascii);
  outFile.close();
}

// [[Rcpp::export]]
void writeOneHaplo(const Rcpp::List& geno,
                   const arma::Col<int>& lociPerChr,
                   arma::uvec lociLoc, int haplo,
                   Rcpp::String filePath, int nThreads){
  arma::Mat<unsigned char> output;
  output = getOneHaplo(geno,lociPerChr,lociLoc,haplo,nThreads);
  std::ofstream outFile;
  outFile.open(filePath, std::ios_base::app);
  output.save(outFile,arma::raw_ascii);
  outFile.close();
}

arma::mat genoToGenoA(const arma::Mat<unsigned char>& geno, 
                      arma::uword ploidy, int nThreads){
  arma::mat output(geno.n_rows,geno.n_cols);
  double dP = double(ploidy);
  arma::vec x(ploidy+1);
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = (double(i)-dP/2.0)*(2.0/dP);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno.n_cols; ++j){
    for(arma::uword i=0; i<geno.n_rows; ++i){
      output(i,j) = x(geno(i,j));
    }
  }
  return output;
}

arma::mat genoToGenoD(const arma::Mat<unsigned char>& geno, 
                      arma::uword ploidy, int nThreads){
  arma::mat output(geno.n_rows,geno.n_cols);
  double dP = double(ploidy);
  arma::vec x(ploidy+1);
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno.n_cols; ++j){
    for(arma::uword i=0; i<geno.n_rows; ++i){
      output(i,j) = x(geno(i,j));
    }
  }
  return output;
}


// Calculates genotype frequency for selected sites
// Intended for use in setEBV with a targetPop
// [[Rcpp::export]]
arma::rowvec calcGenoFreq(const Rcpp::List& geno,
                          const arma::Col<int>& lociPerChr,
                          arma::uvec lociLoc, int nThreads){
  PackedGeno g(geno);
  arma::uvec count = countAlleles(g, lociPerChr, lociLoc, nThreads);
  return freqFromCounts(count, g.nInd, g.ploidy);
}

// Calculates allele frequency on a single chromosome
// [[Rcpp::export]]
arma::rowvec calcChrFreq(const Rcpp::RawVector& geno){
  // A list holding the one chromosome lets it be read in place like any
  // other set of genotypes
  Rcpp::List genoList = Rcpp::List::create(geno);
  PackedGeno g(genoList);
  arma::uword nLoci = g.nBins[0]*8;
  arma::Col<int> lociPerChr(1);
  lociPerChr(0) = nLoci;
  arma::uvec lociLoc(nLoci);
  for(arma::uword i=0; i<nLoci; ++i)
    lociLoc(i) = i+1;
  arma::uvec count = countAlleles(g, lociPerChr, lociLoc, 1);
  return freqFromCounts(count, g.nInd, g.ploidy);
}
