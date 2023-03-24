#include "alphasimr.h"

// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns last element if value is greater than values of x
// Set left to the smallest value of the interval to search
int intervalSearchInt(const arma::Col<int>& x, int value){
  int left = 0;
  
  // Check if crossover is before beginning
  if(x[left]>value){
    return -1;
  }
  int end = x.n_elem-1;
  
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

// Creates the first section of the IBD matrix
arma::Mat<int> firstSection(const arma::Mat<int>& parIbd,
                            int recEnd){
  // Determine last site before new recombination starts
  // Assumes recEnd >= 2 and row 1 of parIbd is always 1
  int matEnd = intervalSearchInt(parIbd.col(1), recEnd-1);
  
  // Return rows for the IBD section
  return parIbd.rows(arma::span(0, matEnd));
}

// Appends to the IBD matrix
arma::Mat<int> middleSection(const arma::Mat<int>& ibd, 
                             const arma::Mat<int>& parIbd,
                             int recStart, 
                             int recEnd){
  // Determine row at or just before recombination
  int matStart = intervalSearchInt(parIbd.col(1), recStart);
  
  // Determine last site before new recombination starts
  // May equal matStart
  int matEnd = intervalSearchInt(parIbd.col(1), recEnd-1);
  
  // Check if previous IBD section matches section after recombination
  if(parIbd(matStart,0) == ibd(ibd.n_rows-1,0)){
    // Increment matStart by 1 to move to next IBD section in parIbd
    // May move past the number of rows in parIbd, handled in next step
    ++matStart;
    
    // Determine if new information is needed from parIbd
    if(matStart>matEnd){
      // No information needed from parIbd
      return ibd;
    }else{
      // Append relevant rows from parIbd to ibd
      return join_cols(ibd, 
                       parIbd.rows(arma::span(matStart, matEnd)));
    }
  }
  
  // Determine if recombination matches start of IBD section (probably rare)
  // Cases where this happens without an IBD change are handled in previous block
  if(parIbd(matStart,1) == recStart){
    // Relevant rows of parIbd can be appended to ibd directly
    return join_cols(ibd, 
                     parIbd.rows(arma::span(matStart, matEnd)));
  }
  
  // Recombination occurred within an IBD block, creating a new segment
  // A new row needs to be added to account for this recombination
  arma::Mat<int> hapStart = { {parIbd(matStart, 0) , recStart} };
  
  // Increment matStart, because first IBD block is now accounted for with hapStart
  ++matStart;
  
  // Determine if additional information is needed from parIbd
  if(matStart>matEnd){
    // No additional information needed from parIbd
    return join_cols(ibd, hapStart);
  }else{
    // Appending additional rows from parIbd to ibd
    return join_cols(join_cols(ibd, hapStart), 
                     parIbd.rows(arma::span(matStart, matEnd)));
  }
}

// Appends to the IBD matrix
arma::Mat<int> lastSection(const arma::Mat<int>& ibd, 
                           const arma::Mat<int>& parIbd,
                           int recStart){
  int matStart = intervalSearchInt(parIbd.col(1), recStart);
  
  // Check if previous IBD section matches section after recombination
  if(parIbd(matStart,0) == ibd(ibd.n_rows-1,0)){
    ++matStart;
    
    // Check if additional information from parIbd is needed
    if(matStart == int(parIbd.n_rows)){
      // No new information needed
      return ibd;
    }else{
      return join_cols(ibd, 
                       parIbd.rows(arma::span(matStart, parIbd.n_rows-1)));
    }
  }
  
  // Check if recombination occurs at new IBD section
  if(parIbd(matStart,1) == recStart){
    // Relevant rows of parIbd can be appended to ibd directly
    return join_cols(ibd, 
                     parIbd.rows(arma::span(matStart, parIbd.n_rows-1)));
  }
  
  // Recombination occurred within an IBD block, creating a new segment
  // A new row needs to be added to account for this recombination
  arma::Mat<int> hapStart = { {parIbd(matStart, 0) , recStart} };
  
  // Increment matStart, because first IBD block is now accounted for with hapStart
  ++matStart;
  
  // Determine if additional information is needed from parIbd
  if(matStart == int(parIbd.n_rows)){
    // No additional information needed from parIbd
    return join_cols(ibd, hapStart);
  }else{
    // Appending additional rows from parIbd to ibd
    return join_cols(join_cols(ibd, hapStart), 
                     parIbd.rows(arma::span(matStart, parIbd.n_rows-1)));
  }
}

// Calculates IBD for individual using recombination data and parental IBD
// [[Rcpp::export]]
arma::field< //chromosome
  arma::field< //ploidy
    arma::Mat<int> > > getNonFounderIbd(
        const arma::field<arma::field<arma::Mat<int> > >& recHist,
        const arma::field<arma::field<arma::Mat<int> > >& mother,
        const arma::field<arma::field<arma::Mat<int> > >& father){
      arma::uword nChr = recHist.n_elem;
      arma::uword motherPloidy = mother(0).n_elem;
      arma::uword fatherPloidy = father(0).n_elem;
      arma::uword ploidy = recHist(0).n_elem;
      arma::field<arma::field<arma::Mat<int> > > output;
      output.set_size(nChr);
      
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
      
      for(arma::uword i=0; i<nChr; ++i){
        
        // Set chromosome to ploidy level
        output(i).set_size(ploidy);
        for(arma::uword j=0; j<ploidy; ++j){
          
          
          arma::uword nCO = recHist(i)(j).n_rows - 1;
          if(j<motherPloidy/2){
            // Pull from mother
            if(nCO==0){
              // Direct copy of chromosome IBD
              output(i)(j) = mother(i)(recHist(i)(j)(0,0)-1);
            }else{
              // Transfer first section
              arma::Mat<int> X = firstSection(mother(i)(recHist(i)(j)(0,0)-1), 
                                              recHist(i)(j)(1,1));
              
              // Transfer middle section(s)?
              if(nCO>1){
                for(arma::uword k=1; k<nCO; ++k){
                  X = middleSection(X, 
                                    mother(i)(recHist(i)(j)(k,0)-1), 
                                    recHist(i)(j)(k,1), 
                                    recHist(i)(j)(k+1,1));
                }
              }
              
              // Transfer last section
              X = lastSection(X,
                              mother(i)(recHist(i)(j)(nCO,0)-1),
                              recHist(i)(j)(nCO,1));
              
              output(i)(j) = X;
            }
          }else{
            // Pull from father
            if(nCO==0){
              // Direct copy of chromosome IBD
              output(i)(j) = father(i)(recHist(i)(j)(0,0)-1 );
            }else{
              // Transfer first section
              arma::Mat<int> X = firstSection(father(i)(recHist(i)(j)(0,0)-1), 
                                              recHist(i)(j)(1,1));
              
              
              // Transfer middle section(s)?
              if(nCO>1){
                for(arma::uword k=1; k<nCO; ++k){
                  X = middleSection(X, 
                                    father(i)(recHist(i)(j)(k,0)-1), 
                                    recHist(i)(j)(k,1), 
                                    recHist(i)(j)(k+1,1));
                }
              }
              
              // Transfer last section
              X = lastSection(X,
                              father(i)(recHist(i)(j)(nCO,0)-1),
                              recHist(i)(j)(nCO,1));
              
              output(i)(j) = X;
            }
          }
        }
      }
      
      return output;
      
    }

// Converts founder IBD into hap format for SimParam
// [[Rcpp::export]]
arma::field< //individual
  arma::field< //chromosome
    arma::field< //ploidy
      arma::Mat<int> > > > getFounderIbd(const arma::field<arma::ivec>& founder, 
                                         arma::uword nChr){
        
        // Generate object for output
        arma::field< //individual
          arma::field< //chromosome
            arma::field< //ploidy
              arma::Mat<int> > > > output;
        output.set_size(founder.n_elem);
        
        // Allocate output object
        for(arma::uword i=0; i<founder.n_elem; ++i){
          output(i).set_size(nChr);
          for(arma::uword j=0; j<nChr; ++j){
            output(i)(j).set_size(founder(i).n_elem);
          }
          for(arma::uword k=0; k<founder(i).n_elem; ++k){
            arma::Mat<int> tmp(1,2,arma::fill::ones);
            tmp(0,0) = founder(i)(k);
            for(arma::uword j=0; j<nChr; ++j){
              output(i)(j)(k) = tmp;
            }
          }
        }
        return output;
      }


// Calculates IBD for individual using recombination data and parental IBD
// [[Rcpp::export]]
arma::Mat<int> createIbdMat(arma::field<arma::field<arma::field<arma::Mat<int> > > >& ibd,
                            arma::uvec chr,
                            arma::uvec nLoci,
                            arma::uword ploidy,
                            arma::uword nThreads){
  // R to C++
  chr -= 1; 
  arma::uword nChr = chr.n_elem;
  arma::uword nInd = ibd.n_elem;
  arma::uword totLoci = accu(nLoci(chr));
  arma::Mat<int> output(totLoci,nInd*ploidy);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<nInd; ++i){
    for(arma::uword j=0; j<ploidy; ++j){
      arma::uword stop,start=0;
      for(arma::uword k=0; k<nChr; ++k){
        arma::uword nSeg = ibd(i)(chr(k))(j).n_rows;
        if(nSeg>1){
          // First segments
          for(arma::uword l=0; l<(nSeg-1); ++l){
            stop = start + ibd(i)(chr(k))(j)(l+1,1) - ibd(i)(chr(k))(j)(l,1) - 1;
            output.col(i*ploidy+j).rows(start,stop).fill(ibd(i)(chr(k))(j)(l,0));
            start = stop + 1;
          }
        }
        // Last segment
        stop = start + nLoci(chr(k)) - ibd(i)(chr(k))(j)(nSeg-1,1);
        output.col(i*ploidy+j).rows(start,stop).fill(ibd(i)(chr(k))(j)(nSeg-1,0));
        start = stop+1;
      }
    }
  }
  
  return output.t();
}
