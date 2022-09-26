#include "alphasimr.h"

/*
 * Genotype data is stored in a field of cubes.
 * The field has length equal to nChr
 * Each cube has dimensions nLoci/8 by ploidy by nInd
 * Output returned with dimensions nInd by nLoci
 */
// [[Rcpp::export]]
arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::Col<int>& lociPerChr,
                                 arma::uvec lociLoc, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  arma::uword ploidy = geno(0).n_cols;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr),arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=0; p<ploidy; ++p){
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(geno(i)(currentByte,p,ind));
          output(ind,loc1) += (unsigned char) workBits[chrLociLoc(0)%8];
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(geno(i)(currentByte,p,ind));
            }
            output(ind,j+loc1) += (unsigned char) workBits[chrLociLoc(j)%8];
          }
        }
      }
    }
  }
  return output;
}

// Extracts genotypes for a specified subset of individuals (indVec)
// The subset does not need to be ordered and is returned in the order submitted
// // [[Rcpp::export]]
// arma::Mat<unsigned char> getGenoSubset(const arma::field<arma::Cube<unsigned char> >& geno, 
//                                        arma::uvec indVec,
//                                        const arma::Col<int>& lociPerChr,
//                                        arma::uvec lociLoc, int nThreads){
//   // R to C++ index correction
//   lociLoc -= 1;
//   
//   arma::uword nInd = indVec.n_elem;
//   arma::uword nChr = geno.n_elem;
//   arma::uword ploidy = geno(0).n_cols;
//   if(nInd < static_cast<arma::uword>(nThreads) ){
//     nThreads = nInd;
//   }
//   arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr),arma::fill::zeros);
//   int loc1;
//   int loc2 = -1;
//   for(arma::uword i=0; i<nChr; ++i){
//     if(lociPerChr(i)>0){
//       // Get loci locations
//       loc1 = loc2+1;
//       loc2 += lociPerChr(i);
//       arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static) num_threads(nThreads)
// #endif
//       for(arma::uword ind=0; ind<nInd; ++ind){
//         std::bitset<8> workBits;
//         arma::uword currentByte, newByte;
//         for(arma::uword p=0; p<ploidy; ++p){
//           currentByte = chrLociLoc(0)/8;
//           workBits = toBits(geno(i)(currentByte,p,indVec(ind)));
//           output(ind,loc1) += (unsigned char) workBits[chrLociLoc(0)%8];
//           for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
//             newByte = chrLociLoc(j)/8;
//             if(newByte != currentByte){
//               currentByte = newByte;
//               workBits = toBits(geno(i)(currentByte,p,indVec(ind)));
//             }
//             output(ind,j+loc1) += (unsigned char) workBits[chrLociLoc(j)%8];
//           }
//         }
//       }
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
arma::Mat<unsigned char> getMaternalGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  arma::uword ploidy = geno(0).n_cols;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr),arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=0; p<ploidy/2; ++p){
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(geno(i)(currentByte,p,ind));
          output(ind,loc1) += (unsigned char) workBits[chrLociLoc(0)%8];
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(geno(i)(currentByte,p,ind));
            }
            output(ind,j+loc1) += (unsigned char) workBits[chrLociLoc(j)%8];
          }
        }
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
arma::Mat<unsigned char> getPaternalGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  arma::uword ploidy = geno(0).n_cols;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr),arma::fill::zeros);
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=ploidy/2; p<ploidy; ++p){
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(geno(i)(currentByte,p,ind));
          output(ind,loc1) += (unsigned char) workBits[chrLociLoc(0)%8];
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(geno(i)(currentByte,p,ind));
            }
            output(ind,j+loc1) += (unsigned char) workBits[chrLociLoc(j)%8];
          }
        }
      }
    }
  }
  return output;
}

// Returns haplotype data in a matrix of nInd*ploidy by nLoci
// [[Rcpp::export]]
arma::Mat<unsigned char> getHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                  const arma::Col<int>& lociPerChr,
                                  arma::uvec lociLoc, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  arma::uword ploidy = geno(0).n_cols;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::Mat<unsigned char> output(nInd*ploidy,arma::sum(lociPerChr));
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=0; p<ploidy; ++p){
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(geno(i)(currentByte,p,ind));
          output(ind*ploidy+p,loc1) = (unsigned char) workBits[chrLociLoc(0)%8];
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(geno(i)(currentByte,p,ind));
            }
            output(ind*ploidy+p,j+loc1) = (unsigned char) workBits[chrLociLoc(j)%8];
          }
        }
      }
    }
  }
  return output;
}

// Returns haplotype data in a matrix of nInd by nLoci for a single
// chromosome group. i.e. just female or male chromosomes for diploids
// [[Rcpp::export]]
arma::Mat<unsigned char> getOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                     const arma::Col<int>& lociPerChr,
                                     arma::uvec lociLoc, int haplo, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  haplo -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  arma::Mat<unsigned char> output(nInd,arma::sum(lociPerChr));
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        currentByte = chrLociLoc(0)/8;
        workBits = toBits(geno(i)(currentByte,haplo,ind));
        output(ind,loc1) = (unsigned char) workBits[chrLociLoc(0)%8];
        for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
          newByte = chrLociLoc(j)/8;
          if(newByte != currentByte){
            currentByte = newByte;
            workBits = toBits(geno(i)(currentByte,haplo,ind));
          }
          output(ind,j+loc1) = (unsigned char) workBits[chrLociLoc(j)%8];
        }
      }
    }
  }
  return output;
}

// Manually sets haplotype data
// [[Rcpp::export]]
arma::field<arma::Cube<unsigned char> > setHaplo(arma::field<arma::Cube<unsigned char> > geno,
                                                 const arma::Mat<unsigned char>& haplo,
                                                 const arma::Col<int>& lociPerChr,
                                                 arma::uvec lociLoc, int nThreads){
  // R to C++ index correction
  lociLoc -= 1;
  
  arma::uword nInd = geno(0).n_slices;
  arma::uword nChr = geno.n_elem;
  arma::uword ploidy = geno(0).n_cols;
  if(nInd < static_cast<arma::uword>(nThreads) ){
    nThreads = nInd;
  }
  int loc1;
  int loc2 = -1;
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){
      // Get loci locations
      loc1 = loc2+1;
      loc2 += lociPerChr(i);
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
      for(arma::uword ind=0; ind<nInd; ++ind){
        std::bitset<8> workBits;
        arma::uword currentByte, newByte;
        for(arma::uword p=0; p<ploidy; ++p){
          currentByte = chrLociLoc(0)/8;
          workBits = toBits(geno(i)(currentByte,p,ind));
          workBits[chrLociLoc(0)%8] = haplo(ind*ploidy+p,loc1);
          geno(i)(currentByte,p,ind) = toByte(workBits);
          for(arma::uword j=1; j<chrLociLoc.n_elem; ++j){
            newByte = chrLociLoc(j)/8;
            if(newByte != currentByte){
              currentByte = newByte;
              workBits = toBits(geno(i)(currentByte,p,ind));
            }
            workBits[chrLociLoc(j)%8] = haplo(ind*ploidy+p,j+loc1);
            geno(i)(currentByte,p,ind) = toByte(workBits);
          }
        }
      }
    }
  }
  return geno;
}

// [[Rcpp::export]]
void writeGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
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
void writeOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
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
arma::rowvec calcGenoFreq(const arma::field<arma::Cube<unsigned char> >& geno, 
                          const arma::Col<int>& lociPerChr,
                          arma::uvec lociLoc, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> tmp = getGeno(geno,
                                         lociPerChr,
                                         lociLoc, nThreads);
  arma::rowvec output = arma::mean(arma::conv_to<arma::mat>::from(tmp),
                                   0)/ploidy;
  return output;
}

// Calculates allele frequency on a single chromsome
// [[Rcpp::export]]
arma::rowvec calcChrFreq(const arma::Cube<unsigned char>& geno){
  arma::field<arma::Cube<unsigned char> > genoList(1);
  genoList(0) = geno;
  arma::uword nLoci = geno.n_rows*8;
  arma::uword ploidy = geno.n_cols;
  arma::Col<int> lociPerChr(1);
  lociPerChr(0) = nLoci;
  arma::uvec lociLoc(nLoci);
  for(arma::uword i=0; i<nLoci; ++i)
    lociLoc(i) = i+1;
  arma::Mat<unsigned char> tmp = getGeno(genoList,
                                         lociPerChr,
                                         lociLoc, 1);
  arma::rowvec output = arma::mean(arma::conv_to<arma::mat>::from(tmp),
                                   0)/ploidy;
  return output;
}
