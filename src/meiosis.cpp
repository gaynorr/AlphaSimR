#include "alphasimr.h"
#include "ibd.h"
#include <algorithm>
#include <cstring>
#include <vector>

namespace {

// Samples the position of the first chiasma at or after position 0 for a
// stationary gamma renewal process, and optionally the position of the
// chiasma preceding position 0.
//
// A stationary renewal process is not started by burning in from a position
// several Morgans upstream. The renewal interval that covers position 0 is
// length biased, which for Gamma(shape, scale) spacings is exactly
// Gamma(shape+1, scale), and position 0 falls uniformly within that
// interval. Splitting the interval at a uniform point therefore gives the
// backward and forward recurrence times directly, using two random deviates
// instead of the twenty or so consumed by a burn-in.
//
// backward, when not null, receives the position of the preceding chiasma.
// It is zero or negative.
double sampleFirstChiasma(double shape, double scale,
                          alphasimrRng::rngEngine& rng,
                          double* backward=nullptr){
  double interval = alphasimrRng::gammaSampler(shape+1.0, scale)(rng);
  double u = alphasimrRng::runif(rng);
  if(backward != nullptr){
    *backward = -u*interval;
  }
  return (1.0-u)*interval;
}

// Creates one stable dqrng substream for each work item of the current call.
//
// Streams are built by walking: each one is a single long jump on from the
// last. dqrng applies clone(stream) as a loop of stream long jumps, so
// asking for stream k costs k jumps and building n streams by index would
// cost O(n^2). Walking costs one jump per stream and reaches the same
// states. A long jump is 2^192 draws, far more than a work item consumes,
// so the streams do not overlap.
//
// Stream ids follow work item order, which is fixed by the chromosome and
// block indices, so results stay reproducible if the OpenMP thread count
// changes.
std::vector<alphasimrRng::rngPtr> makeWorkRngs(arma::uword nWork) {
  std::vector<alphasimrRng::rngPtr> workRngs;
  if (nWork == 0) {
    return workRngs;
  }
  dqrng::rng64_t baseRng = alphasimrRng::createRng();
  workRngs.reserve(nWork);
  workRngs.push_back(baseRng->clone(1));
  for (arma::uword i = 1; i < nWork; ++i) {
    workRngs.push_back(workRngs[i - 1]->clone(1));
  }
  return workRngs;
}

} // namespace

// Class for storing recombination history
class RecHist{
public:
  arma::field< //individual
    arma::field< //chromosome
      arma::field< //ploidy
        arma::Mat<int> > > > hist; //(chr, site)
  
  // Allocates space for history
  void setSize(arma::uword nInd, 
               arma::uword nChr, 
               arma::uword ploidy);
  
  // Append new recombination history
  void addHist(arma::Mat<int>& input, 
               arma::uword nInd, 
               arma::uword chrGroup,
               arma::uword chrInd);
  
  // Internal accessor with no current callers or R interface;
  // R code reads histories through pack(),
  // so R-level tests cannot reach getHist().
  arma::Mat<int> getHist(arma::uword ind, 
                         arma::uword chr,
                         arma::uword par);

  // Packs each individual's history into one integer vector for R
  Rcpp::List pack();
};

// Allocates space for history
void RecHist::setSize(arma::uword nInd, 
                      arma::uword nChr, 
                      arma::uword ploidy=2){
  hist.set_size(nInd);
  for(arma::uword i=0; i<nInd; ++i){
    hist(i).set_size(nChr);
    for(arma::uword j=0; j<nChr; ++j){
      hist(i)(j).set_size(ploidy);
    }
  }
}

// Append new recombination history
void RecHist::addHist(arma::Mat<int>& input, 
             arma::uword nInd, 
             arma::uword chrGroup,
             arma::uword chrInd){
  hist(nInd)(chrGroup)(chrInd) = input;
}

// Access history
arma::Mat<int> RecHist::getHist(arma::uword ind, 
                                arma::uword chr,
                                arma::uword par){
  return hist(ind)(chr)(par);
}

// Packs each individual's history into one integer vector for R. Returning
// the nested field made R allocate about thirty objects per individual.
Rcpp::List RecHist::pack(){
  Rcpp::List output(hist.n_elem);
  for(arma::uword i=0; i<hist.n_elem; ++i){
    output[i] = packRecord(hist(i));
  }
  return output;
}

// Samples the locations for chiasmata via a gamma process
// end, the length of the interval used to sample
// v, the interference parameter
// p, the proportion of non-interfering crossovers
// rng, the explicit dqrng stream used for all random draws in this call
//
// The type 1 (interfering) chiasmata come from a stationary gamma renewal
// process. The first one is drawn from the equilibrium distribution by
// sampleFirstChiasma, so only the deviates that fall on the chromosome are
// sampled. The type 2 (non-interfering) chiasmata come from a Poisson
// process, which is memoryless and so needs no equilibrium start.
arma::vec sampleChiasmata(double end, double v,
                          double p, alphasimrRng::rngEngine& rng){
  if((1.0-p)<1.0e-6){
    // No crossover interference
    // Switching to count-location model
    arma::uword n = alphasimrRng::samplePoisson(2.0*end, rng);
    // Positions are uniform over the chromosome, so the deviates have to be
    // scaled by its genetic length. Without the scale every chiasma lands in
    // the first Morgan.
    arma::vec x = alphasimrRng::runifVec(n, rng)*end;
    return sort(x);

  }else{
    // Using gamma or gamma-sprinkling model
    std::vector<double> pos;
    pos.reserve(static_cast<std::size_t>(2.0*end)+8);

    // Sample type 1 chiasmata
    double scale;
    if(p<1.0e-6){ // Gamma model
      scale = 1.0/(2.0*v);
    }else{ // Gamma sprinkling model
      scale = 1.0/(2.0*v*(1.0-p));
    }
    alphasimrRng::gammaSampler sampleGap(v, scale);
    double x = sampleFirstChiasma(v, scale, rng);
    while(x<end){
      pos.push_back(x);
      x += sampleGap(rng);
    }

    if(!(p<1.0e-6)){ // Gamma sprinkling model
      // Sample type 2 chiasmata from a Poisson process
      alphasimrRng::gammaSampler sampleSprinkle(1.0, 1.0/(2.0*p));
      x = sampleSprinkle(rng);
      while(x<end){
        pos.push_back(x);
        x += sampleSprinkle(rng);
      }

      // Combine type 1 and type 2 crossovers in order
      std::sort(pos.begin(), pos.end());
    }

    if(pos.empty()){
      return arma::vec();
    }
    return arma::vec(pos.data(), pos.size());
  }
}

// Samples the locations for chiasmata via a gamma process for a quadrivalent
// CO interference is assumed to occur between all arms
// The first arm is sampled at random
// exchange, the positions where chromosomes switch
// end, the length of the interval used to sample
// v, the interference parameter
// p, the proportion of non-interfering crossovers
// rng, the explicit dqrng stream used for all random draws in this call
// n2, the number of gamma deviates sampled for all other arms
arma::field<arma::vec> sampleQuadChiasmata(double exchange, double end, double v,
                                           double p, alphasimrRng::rngEngine& rng,
                                           arma::uword n2=8){
  arma::field<arma::vec> output(4);

  // Randomly set order of chromosome arms
  arma::uvec arm = {0, 1, 2, 3};
  alphasimrRng::shuffle(arm, rng);
  double nearest, terminator, prob;
  
  if((1.0-p)<1.0e-6){
    // All chiasmata from type 2 pathway
    // Changing v and p to model type 2 with type 1 pathway
    p = 0.0;
    v = 1.0;
  }
  
  // First arm
  // Sampled from a stationary gamma renewal process. The chiasma preceding
  // position 0 is retained so that nearest still measures the true distance
  // back to the last chiasma when the arm itself contains none.
  if(arm(0)%2){ // Tail
    terminator = end - exchange;
  }else{ // Head
    terminator = exchange;
  }
  {
    double scale = 1.0/(2.0*v*(1.0-p));
    alphasimrRng::gammaSampler sampleGap(v, scale);
    double previous;
    double x = sampleFirstChiasma(v, scale, rng, &previous);
    std::vector<double> pos;
    pos.reserve(static_cast<std::size_t>(2.0*terminator)+8);
    pos.push_back(previous);
    while(x<terminator){
      pos.push_back(x);
      x += sampleGap(rng);
    }
    output(arm(0)) = arma::vec(pos.data(), pos.size());
  }
  nearest = terminator - output(arm(0))(output(arm(0)).n_elem-1);
  output(arm(0)) = output(arm(0))(find(output(arm(0))>0));
  if(arm(0)%2){ // Tail
    output(arm(0)) = sort(end-output(arm(0)));
  }
  
  // All other arms
  for(arma::uword i=1; i<4; ++i){
    output(arm(i)).set_size(1+n2);
    prob = R::pgamma(nearest, v, 1.0/(2.0*v*(1.0-p)), 1, 0);
    double u = alphasimrRng::runif(rng);
    u = u*(1.0-prob)+prob;
    output(arm(i))(0) = R::qgamma(u, v, 1.0/(2.0*v*(1.0-p)), 1, 0) - nearest;
    if(output(arm(i))(0) < nearest){
      nearest = output(arm(i))(0);
    }
    output(arm(i))(arma::span(1,n2)) = alphasimrRng::rgammaVec(
      n2, v, 1.0/(2.0*v*(1.0-p)), rng);
    output(arm(i)) = cumsum(output(arm(i)));
    if(arm(i)%2){ // Tail
      terminator = end - exchange;
    }else{ // Head
      terminator = exchange;
    }
    while( output(arm(i))(output(arm(i)).n_elem-1) < terminator ){
      arma::vec tmp = alphasimrRng::rgammaVec(n2, v, 1.0/(2.0*v*(1.0-p)), rng);
      tmp = cumsum(tmp) + output(arm(i))(output(arm(i)).n_elem-1);
      output(arm(i)) = join_cols(output(arm(i)), tmp);
    }
    output(arm(i)) = output(arm(i))(find(output(arm(i))<terminator));
    if(arm(i)%2){ // Tail
      output(arm(i)) += exchange;
    }else{ // Head
      output(arm(i)) = sort(exchange-output(arm(i)));
    }
  }
  
  if(p>1.0e-6){ // Sprinkle recombinations
    for(arma::uword i=0; i<4; ++i){
      if(arm(i)%2){ // Tail
        terminator = end - exchange;
      }else{ // Head
        terminator = exchange;
      }
      
      // Sample type 2 deviates from a gamma distribution
      arma::vec type2 = alphasimrRng::rgammaVec(n2, 1.0, 1.0/(2.0*p), rng);
      
      // Find locations on genetic map
      type2 = cumsum(type2);
      
      // Add additional values if max position less than terminator
      while(type2(type2.n_elem-1)<terminator){
        arma::vec tmp = alphasimrRng::rgammaVec(n2, 1.0, 1.0/(2.0*p), rng);
        tmp = cumsum(tmp) + type2(type2.n_elem-1);
        type2 = join_cols(type2, tmp);
      }
      
      // Select values less than the end
      type2 = type2(find(type2<terminator));
      
      // Add exchange to genetic position if in tail
      if(arm(i)%2){ // Tail
        type2 += exchange;
      }
      
      // Combine type 1 and 2 crossovers and sort
      output(arm(i)) = sort(join_cols(output(arm(i)), type2));
    }
  }
  
  return output;
}


// Counts the loci in genMap whose positions are at or before value. The
// search begins at start, which must not exceed the count, so a caller
// resolving crossovers in increasing order can pass the previous result.
//
// Genetic map positions are absolute, so a crossover can fall before the
// first locus, giving a count of zero, or after the last, giving the number
// of loci. Neither is an error. The count plus one is the 1-indexed site at
// which the segment following the crossover begins, which is the form a
// recombination history records.
arma::uword countLociBefore(const arma::vec& genMap, double value,
                            arma::uword start=0){
  return static_cast<arma::uword>(
    std::upper_bound(genMap.begin()+start, genMap.end(), value) -
      genMap.begin());
}

// Removes hidden crossovers from recombination map
// Assumes the first row is site 1 and that sites never decrease. Several
// rows can share site 1, because every crossover before the first locus is
// recorded there, and the last of them sets the starting chromosome.
arma::Mat<int> removeDoubleCO(const arma::Mat<int>& X){
  if(X.n_rows<2){
    return X;
  }

  // Initially assume all rows are useful
  arma::Col<int> take(X.n_rows,arma::fill::ones);

  // Remove unobserved crossovers (site doesn't change)
  // Works backward, because the last crossover is observed. Row 0 is
  // included so that crossovers before the first locus collapse into it.
  for(arma::uword i=(X.n_rows-1); i>0; --i){
    if(X(i-1,1) == X(i,1)){
      take(i-1) = 0;
    }
  }
  
  // Remove redundant records (chromosome doesn't change)
  // Starts from the first row kept, which is the last row at site 1
  arma::uword first = 0;
  while(take(first) == 0){
    ++first;
  }
  int lastChr = X(first,0);
  for(arma::uword i=first+1; i<X.n_rows; ++i){
    if(take(i) == 1){
      if(X(i,0) == lastChr){
        take(i) = 0;
      }else{
        lastChr = X(i,0);
      }
    }
  }
  return X.rows(find(take>0));
}

// Finds recombination map for a bivalent pair
// genMap, chromosome genetic map
// genLen, genetic length of the chromosome
// v, the interference parameter
// p, the proportion of non-interfering crossovers
// rng, the explicit dqrng stream used for all random draws in this call
//
// Map positions are absolute, so the chromosome runs from position 0 to
// genLen, which can lie past the last locus, and chiasmata are sampled over
// all of it. Those falling before the first locus are resolved to site 1,
// and those after the last locus affect no locus and are not recorded.
arma::Mat<int> findBivalentCO(const arma::vec& genMap, double genLen,
                              double v, double p,
                              alphasimrRng::rngEngine& rng){
  arma::uword nBefore=0, readChr=0, nCO;
  
  // Find crossover positions
  arma::vec posCO = sampleChiasmata(genLen, v, p, rng);
  if(posCO.n_elem==0){
    arma::Mat<int> output(1,2,arma::fill::ones);
    return output;
  }
  
  // Thin crossovers 
  arma::vec thin = alphasimrRng::runifVec(posCO.n_elem, rng);
  posCO = posCO(find(thin>0.5));
  nCO = posCO.n_elem;
  
  arma::Mat<int> output(nCO+1,2);
  if(nCO==0){
    output.ones();
    return output;
  }
  
  // Find crossover sites on map
  // Crossovers are sorted, so once one falls after the last locus every
  // later one does too
  output.row(0).ones();
  arma::uword nRows = 1;
  for(arma::uword i=0; i<nCO; ++i){
    ++readChr;
    readChr = readChr%2;
    nBefore = countLociBefore(genMap,posCO(i),nBefore);
    if(nBefore==genMap.n_elem){
      break;
    }
    output(nRows,0) = readChr+1;
    output(nRows,1) = nBefore+1;
    ++nRows;
  }
  
  return removeDoubleCO(output.head_rows(nRows));
}

/*
 * Finds recombination maps for a quadrivalent "cross-type" configuration
 * The configuration for chromosome pairing is as follows:
 *  Arm 0: chromosome heads 1 and 2
 *  Arm 1: chromosome tail 2 and 3
 *  Arm 2: chromosome heads 3 and 4
 *  Arm 3: chromosome tails 1 and 4
 * The exchange point between pairings is sampled at random
 * A centromere from the first chromosome is always selected
 * The second centromere is sampled at random
 * rng is the explicit dqrng stream used for all random draws in this call
 */
arma::field<arma::Mat<int> > findQuadrivalentCO(const arma::vec& genMap,
                                                double genLen,
                                                double centromere, double v,
                                                double p,
                                                alphasimrRng::rngEngine& rng){
  arma::field<arma::Mat<int> > output(2);
  
  // Sample the exchange point
  double exchange = alphasimrRng::runif(rng) * genLen;
  
  // Determine crossover positions
  // Returns field with crossover positions in each arm of the quadrivalent
  arma::field<arma::vec> posCO = sampleQuadChiasmata(exchange, genLen, v, p, rng);
  
  // Set chromatid configuration for each chiasmata
  arma::field<arma::umat> chromatidPairs(4); // matches posCO
  for(arma::uword i=0; i<4; ++i){
    // Create table for chromatid pairs on an arm
    chromatidPairs(i).set_size(posCO(i).n_elem,2);
    
    // Assign pairs if there are chiasmata
    if(chromatidPairs(i).n_rows>0){
      // Initializing with "0" chromatid
      chromatidPairs(i).zeros();
      
      // Randomly switch to "1" chromatid
      for(arma::uword j=0; j<chromatidPairs(i).n_elem; ++j){
        if(alphasimrRng::runif(rng)>0.5){
          chromatidPairs(i).at(j) = 1;
        }
      }
    }
    
  }
  
  // Allocate output with a naive maximum number of COs
  arma::uword maxCO=0;
  for(arma::uword i=0; i<4; ++i){
    maxCO = std::max(maxCO, posCO(i).n_elem);
  }
  maxCO *= 2;
  output(0).set_size(maxCO+1,2);
  output(1).set_size(maxCO+1,2);
  
  // Select centromeres (which chromosome and chromatid)
  // Always taking chromosome 1 (1-4) and chromatid 1 (0-1) 
  arma::uvec chromosome(2, arma::fill::ones);
  arma::uvec chromatid(2, arma::fill::ones);
  chromosome(1) = alphasimrRng::sampleInt(1,3, rng)(0) + 2; // 2-4
  chromatid(1) = alphasimrRng::sampleInt(1,2, rng)(0); // 0-1
  
  // Find starting chromosomes and chromatids by working backward
  // from selected centromeres to start of chromosome (head)
  arma::uword currentChromosome, currentChromatid;
  for(arma::uword i=0; i<2; ++i){
    // Identify starting chromosome and chromatid
    currentChromosome = chromosome(i);
    currentChromatid = chromatid(i);
    
    if(exchange>centromere){ // Centromere is in the head 
      
      if(currentChromosome<3){ // currentChromosome is 1 or 2
        
        // Loop through all chiasmata on arm 0
        for(arma::uword j=posCO(0).n_elem; j>0; --j){
          
          // Check if chiasmata is before centromere, ignore if not
          if(posCO(0)(j-1)<centromere){
            
            switch(currentChromosome){
            case 1: // Check if there's a switch between chr 1 and 2
              if(chromatidPairs(0)(j-1,0) == currentChromatid){
                currentChromosome = 2;
                currentChromatid = chromatidPairs(0)(j-1,1);
              }
              break;
            case 2: // Check if there is a switch between chr 2 and 1
              if(chromatidPairs(0)(j-1,1) == currentChromatid){
                currentChromosome = 1;
                currentChromatid = chromatidPairs(0)(j-1,0);
              }
            }
          }
        }
        
      }else{ // currentChromosome is 3 or 4
        
        // Loop through all chiasmata on arm 2
        for(arma::uword j=posCO(2).n_elem; j>0; --j){
          
          // Check if chiasmata is before centromere, ignore if not
          if(posCO(2)(j-1)<centromere){
            
            switch(currentChromosome){
            case 3: // Check if there's a switch between chr 3 and 4
              if(chromatidPairs(2)(j-1,0) == currentChromatid){
                currentChromosome = 4;
                currentChromatid = chromatidPairs(2)(j-1,1);
              }
              break;
            case 4: // Check if there's a switch between chr 4 and 3
              if(chromatidPairs(2)(j-1,1) == currentChromatid){
                currentChromosome = 3;
                currentChromatid = chromatidPairs(2)(j-1,0);
              }
            }
          }
        }
      }
      
    }else{ // Centromere is in the tail
      
      if((currentChromosome==1) | (currentChromosome==4)){ // Working on arm 3
        
        // Find chromosome and chromatid before transition by looping 
        // through chiasmata on arm 3
        for(arma::uword j=posCO(3).n_elem; j>0; --j){
          
          // Check if chiasmata is before centromere, ignore if not
          if(posCO(3)(j-1)<centromere){
            switch(currentChromosome){
            case 1: // Check if there's a switch between chr 1 and 4
              if(chromatidPairs(3)(j-1,0) == currentChromatid){
                currentChromosome = 4;
                currentChromatid = chromatidPairs(3)(j-1,1);
              }
              break;
            case 4: // Check if there's a switch between chr 4 and 1
              if(chromatidPairs(3)(j-1,1) == currentChromatid){
                currentChromosome = 1;
                currentChromatid = chromatidPairs(3)(j-1,0);
              }
            }
          }
        }
        
        // Find starting chromosome and chromatid by working back through head
        switch(currentChromosome){
        case 1: // Work through arm 0
          for(arma::uword j=posCO(0).n_elem; j>0; --j){
            switch(currentChromosome){
            case 1: // Check if there's a switch between chr 1 and 2
              if(chromatidPairs(0)(j-1,0) == currentChromatid){
                currentChromosome = 2;
                currentChromatid = chromatidPairs(0)(j-1,1);
              }
              break;
            case 2: // Check if there's a switch between chr 2 and 1
              if(chromatidPairs(0)(j-1,1) == currentChromatid){
                currentChromosome = 1;
                currentChromatid = chromatidPairs(0)(j-1,0);
              }
            }
          }
          break;
        case 4: // Work through arm 2
          for(arma::uword j=posCO(2).n_elem; j>0; --j){
            switch(currentChromosome){
            case 3: // Check if there's a switch between chr 3 and 4
              if(chromatidPairs(2)(j-1,0) == currentChromatid){
                currentChromosome = 4;
                currentChromatid = chromatidPairs(2)(j-1,1);
              }
              break;
            case 4: // Check if there's a switch between chr 4 and 3
              if(chromatidPairs(2)(j-1,1) == currentChromatid){
                currentChromosome = 3;
                currentChromatid = chromatidPairs(2)(j-1,0);
              }
            }
          }
        }
        
      }else{ // Working on arm 1
        
        // Find chromosome and chromatid before transition by looping 
        // through chiasmata on arm 1
        for(arma::uword j=posCO(1).n_elem; j>0; --j){
          
          // Check if chiasmata is before centromere, ignore if not
          if(posCO(1)(j-1)<centromere){
            switch(currentChromosome){
            case 2: // Check if there's a switch between chr 2 and 3
              if(chromatidPairs(1)(j-1,0) == currentChromatid){
                currentChromosome = 3;
                currentChromatid = chromatidPairs(1)(j-1,1);
              }
              break;
            case 3: // Check if there's a switch between chr 3 and 2
              if(chromatidPairs(1)(j-1,1) == currentChromatid){
                currentChromosome = 2;
                currentChromatid = chromatidPairs(1)(j-1,0);
              }
            }
          }
        }
        
        // Find starting chromosome and chromatid by working back through head
        switch(currentChromosome){
        case 2: // Work through arm 0
          for(arma::uword j=posCO(0).n_elem; j>0; --j){
            switch(currentChromosome){
            case 1: // Check if there's a switch between chr 1 and 2
              if(chromatidPairs(0)(j-1,0) == currentChromatid){
                currentChromosome = 2;
                currentChromatid = chromatidPairs(0)(j-1,1);
              }
              break;
            case 2: // Check if there's a switch between chr 2 and 1
              if(chromatidPairs(0)(j-1,1) == currentChromatid){
                currentChromosome = 1;
                currentChromatid = chromatidPairs(0)(j-1,0);
              }
            }
          }
          break;
        case 3: // Work through arm 2
          for(arma::uword j=posCO(2).n_elem; j>0; --j){
            switch(currentChromosome){
            case 3: // Check if there's a switch between chr 3 and 4
              if(chromatidPairs(2)(j-1,0) == currentChromatid){
                currentChromosome = 4;
                currentChromatid = chromatidPairs(2)(j-1,1);
              }
              break;
            case 4: // Check if there's a switch between chr 4 and 3
              if(chromatidPairs(2)(j-1,1) == currentChromatid){
                currentChromosome = 3;
                currentChromatid = chromatidPairs(2)(j-1,0);
              }
            }
          }
        }
      }
    }
    
    // Fill in crossover map by working from head to tail
    arma::uword nBefore=0, nCO=0;
    output(i)(0,0) = currentChromosome;
    output(i)(0,1) = 1;
    
    // Records a crossover onto the chromosome now being read. The arms are
    // read from head to tail, so crossovers arrive in increasing order. One
    // past the last locus affects no locus and is not recorded, and neither
    // is any after it, which keeps every site in the history within the
    // chromosome's loci.
    auto recordCO = [&](double pos){
      nBefore = countLociBefore(genMap,pos,nBefore);
      if(nBefore<genMap.n_elem){
        ++nCO;
        output(i)(nCO,0) = currentChromosome;
        output(i)(nCO,1) = nBefore+1;
      }
    };
    
    if(currentChromosome<3){ // Start in arm 0
      
      // Fill crossovers in the head
      for(arma::uword j=0; j<posCO(0).n_elem; ++j){
        
        // Check if chiasmata involves current chromatid
        switch(currentChromosome){
        case 1: 
          if(chromatidPairs(0)(j,0) == currentChromatid){
            currentChromosome = 2;
            currentChromatid = chromatidPairs(0)(j,1);
            recordCO(posCO(0)(j));
          }
          break;
        case 2:
          if(chromatidPairs(0)(j,1) == currentChromatid){
            currentChromosome = 1;
            currentChromatid = chromatidPairs(0)(j,0);
            recordCO(posCO(0)(j));
          }
        }
      }
      
      // Fill crossovers in the tail
      if(currentChromosome==1){ // Move to arm 3
        for(arma::uword j=0; j<posCO(3).n_elem; ++j){
          switch(currentChromosome){
          case 1:
            if(chromatidPairs(3)(j,0) == currentChromatid){
              currentChromosome = 4;
              currentChromatid = chromatidPairs(3)(j,1);
              recordCO(posCO(3)(j));
            }
            break;
          case 4:
            if(chromatidPairs(3)(j,1) == currentChromatid){
              currentChromosome = 1;
              currentChromatid = chromatidPairs(3)(j,0);
              recordCO(posCO(3)(j));
            }
          }
        }
      }else{ // currentChromosome = 2, move to arm 1
        for(arma::uword j=0; j<posCO(1).n_elem; ++j){
          switch(currentChromosome){
          case 2:
            if(chromatidPairs(1)(j,0) == currentChromatid){
              currentChromosome = 3;
              currentChromatid = chromatidPairs(1)(j,1);
              recordCO(posCO(1)(j));
            }
            break;
          case 3:
            if(chromatidPairs(1)(j,1) == currentChromatid){
              currentChromosome = 2;
              currentChromatid = chromatidPairs(1)(j,0);
              recordCO(posCO(1)(j));
            }
          }
        }
      }
    }else{ // currentChromosome>2, start in arm 2
      
      // Fill crossovers in the head
      for(arma::uword j=0; j<posCO(2).n_elem; ++j){
        switch(currentChromosome){
        case 3:
          if(chromatidPairs(2)(j,0) == currentChromatid){
            currentChromosome = 4;
            currentChromatid = chromatidPairs(2)(j,1);
            recordCO(posCO(2)(j));
          }
          break;
        case 4:
          if(chromatidPairs(2)(j,1) == currentChromatid){
            currentChromosome = 3;
            currentChromatid = chromatidPairs(2)(j,0);
            recordCO(posCO(2)(j));
          }
        }
      }
      
      // Fill crossovers in the tail
      if(currentChromosome==4){ // Move to arm 3
        for(arma::uword j=0; j<posCO(3).n_elem; ++j){
          switch(currentChromosome){
          case 1:
            if(chromatidPairs(3)(j,0) == currentChromatid){
              currentChromosome = 4;
              currentChromatid = chromatidPairs(3)(j,1);
              recordCO(posCO(3)(j));
            }
            break;
          case 4:
            if(chromatidPairs(3)(j,1) == currentChromatid){
              currentChromosome = 1;
              currentChromatid = chromatidPairs(3)(j,0);
              recordCO(posCO(3)(j));
            }
          }
        }
      }else{ // currentChromosome = 3, move to arm 1
        for(arma::uword j=0; j<posCO(1).n_elem; ++j){
          switch(currentChromosome){
          case 2:
            if(chromatidPairs(1)(j,0) == currentChromatid){
              currentChromosome = 3;
              currentChromatid = chromatidPairs(1)(j,1);
              recordCO(posCO(1)(j));
            }
            break;
          case 3:
            if(chromatidPairs(1)(j,1) == currentChromatid){
              currentChromosome = 2;
              currentChromatid = chromatidPairs(1)(j,0);
              recordCO(posCO(1)(j));
            }
          }
        }
      }
    }
    output(i) = output(i).rows(arma::span(0,nCO));
    output(i) = removeDoubleCO(output(i));
  }
  
  return output;
}

// Copies packed genotypes for sites [start, stop) from inChr to outChr.
// start and stop are 1-indexed site numbers, as recorded in a recombination
// history. Both haplotypes hold nBytes bytes.
//
// The haplotypes are passed as pointers rather than as arma::Col so that a
// caller can hand over a column of a cube without Armadillo materializing a
// temporary copy of it.
void transferGeno(const unsigned char* inChr,
                  unsigned char* outChr,
                  arma::uword nBytes,
                  int start,
                  int stop){
  start -= 1; // R to C++
  stop -= 1; // R to C++
  std::bitset<8> inBits, outBits;
  int startByte = start / 8;
  int stopByte = stop / 8;
  int startBit = start % 8;
  int stopBit = stop % 8;
  // Transfer partial start
  if(startBit != 0){
    inBits = toBits(inChr[startByte]);
    outBits = toBits(outChr[startByte]);
    if(stopByte > startByte){
      // Transferring more than this byte
      for(int i=startBit; i<8; ++i){
        outBits[i] = inBits[i];
      }
      outChr[startByte] = toByte(outBits);
      startBit = 0;
      ++startByte;
    }else{
      // Only transferring within this byte
      for(int i=startBit; i<stopBit; ++i){
        outBits[i] = inBits[i];
      }
      outChr[startByte] = toByte(outBits);
      return;
    }
  }
  // Transfer full bytes
  if(stopByte >  startByte){
    std::memcpy(outChr+startByte, inChr+startByte,
                static_cast<std::size_t>(stopByte-startByte));
    startByte = stopByte;
  }
  // Transfer partial stop
  if(nBytes == static_cast<arma::uword>(startByte) ){
    // End has been reached
    return;
  }else{
    if(stopBit > startBit){
      inBits = toBits(inChr[startByte]);
      outBits = toBits(outChr[startByte]);
      for(int i = startBit; i<stopBit; ++i){
        outBits[i] = inBits[i];
      }
      outChr[startByte] = toByte(outBits);
    }
  }
}

// Simulates a gamete using a count-location model for recombination
// rng is the explicit dqrng stream used for crossover sampling and thinning.
//
// The parental haplotypes and the gamete are passed as pointers to nBins
// bytes of packed genotypes. Passing pointers lets the caller read parents
// straight out of a genotype cube and write the gamete straight into the
// progeny cube, with no intermediate copies of either.
void bivalent(const unsigned char* chr1,
              const unsigned char* chr2,
              arma::uword nBins,
              const arma::vec& genMap,
              double genLen,
              double v,
              double p,
              unsigned char* output,
              arma::Mat<int>& hist,
              alphasimrRng::rngEngine& rng){
  hist = findBivalentCO(genMap, genLen, v, p, rng);
  if(hist.n_rows==1){
    // No crossovers between loci, so the gamete is a copy of one chromosome.
    // It is not always chromosome 1, because crossovers before the first
    // locus change which chromosome the gamete starts on.
    std::memcpy(output, (hist(0,0)==1) ? chr1 : chr2, nBins);
  }else{
    int nSites = int(nBins)*8+1;

    // Fill-in based on recombination history
    for(arma::uword i=0; i<(hist.n_rows-1); ++i){
      switch(hist(i,0)){
      case 1: //Chromosome 1
        transferGeno(chr1, output, nBins,
                     hist(i,1), hist(i+1,1));
        break;
      case 2: //Chromosome 2
        transferGeno(chr2, output, nBins,
                     hist(i,1), hist(i+1,1));
      }
    }

    // Fill-in last sites
    switch(hist(hist.n_rows-1,0)){
    case 1:
      transferGeno(chr1, output, nBins,
                   hist(hist.n_rows-1,1),
                   nSites);
      break;
    case 2:
      transferGeno(chr2, output, nBins,
                   hist(hist.n_rows-1,1),
                   nSites);
    }
  }
}

// Resolves one gamete of a quadrivalent from its recombination history.
// chr holds pointers to the four parental haplotypes, each of nBins bytes,
// indexed by the chromosome numbers (1-4) recorded in the history.
void resolveQuadGamete(const unsigned char* const chr[4],
                       arma::uword nBins,
                       const arma::Mat<int>& hist,
                       unsigned char* output){
  if(hist.n_rows==1){
    // No crossovers, so the gamete is a copy of a single chromosome
    std::memcpy(output, chr[hist(0,0)-1], nBins);
    return;
  }

  // Fill-in based on recombination history
  for(arma::uword i=0; i<(hist.n_rows-1); ++i){
    transferGeno(chr[hist(i,0)-1], output, nBins,
                 hist(i,1), hist(i+1,1));
  }

  // Fill-in last sites
  transferGeno(chr[hist(hist.n_rows-1,0)-1], output, nBins,
               hist(hist.n_rows-1,1),
               int(nBins)*8+1);
}

// Simulates a pair of gametes using a count-location model for recombination
// rng is the explicit dqrng stream used for crossover sampling and thinning.
//
// As in bivalent, the parental haplotypes and the gametes are passed as
// pointers to nBins bytes of packed genotypes so that no copies are made on
// the way in or out.
void quadrivalent(const unsigned char* chr1,
                  const unsigned char* chr2,
                  const unsigned char* chr3,
                  const unsigned char* chr4,
                  arma::uword nBins,
                  const arma::vec& genMap,
                  double genLen,
                  double centromere,
                  double v,
                  double p,
                  unsigned char* output1,
                  unsigned char* output2,
                  arma::Mat<int>& hist1,
                  arma::Mat<int>& hist2,
                  alphasimrRng::rngEngine& rng){
  const unsigned char* const chr[4] = {chr1, chr2, chr3, chr4};

  arma::field<arma::Mat<int> > output;
  output = findQuadrivalentCO(genMap, genLen, centromere, v, p, rng);

  hist1 = output(0);
  hist2 = output(1);

  resolveQuadGamete(chr, nBins, hist1, output1);
  resolveQuadGamete(chr, nBins, hist2, output2);
}

// Makes crosses between diploid individuals.
// motherGeno: female genotypes
// mother: female parents
// fatherGeno: male genotypes
// father: male parents
// femaleMap: chromosome genetic maps
// maleMap: chromosome genetic maps
// femaleGenLen: chromosome genetic lengths for the female map
// maleGenLen: chromosome genetic lengths for the male map
// trackRec: track recombination
// motherPloidy: ploidy level of mother 
// fatherPloidy: ploidy level of father
// v: interference parameter for gamma model
// p: proportion of non-interfering crossovers
// quadProb: probability of quadrivalent formation
// nThreads: number of threads for parallel computing
// [[Rcpp::export]]
Rcpp::List cross(
    const arma::field<arma::Cube<unsigned char> >& motherGeno, 
    arma::uvec mother,
    const arma::field<arma::Cube<unsigned char> >& fatherGeno, 
    arma::uvec father,
    const arma::field<arma::vec>& femaleMap,
    const arma::field<arma::vec>& maleMap,
    const arma::vec& femaleGenLen,
    const arma::vec& maleGenLen,
    bool trackRec,
    arma::uword motherPloidy,
    arma::uword fatherPloidy,
    double v,
    double p,
    const arma::vec& motherCentromere,
    const arma::vec& fatherCentromere,
    double quadProb,
    int nThreads){
  // The parents pick the genotype slices that meiosis reads through raw
  // pointers, so they are checked here, before the parallel region. The
  // indexes are still one based, and a zero or negative index from R has
  // become 0 or a very large value.
  if(mother.n_elem!=father.n_elem){
    Rcpp::stop("mother and father must have the same length");
  }
  if((motherGeno.n_elem==0) || (fatherGeno.n_elem!=motherGeno.n_elem)){
    Rcpp::stop("mother and father genotypes must have the same chromosomes");
  }
  // Each chromosome's genetic length is read inside the parallel region
  if((femaleGenLen.n_elem!=motherGeno.n_elem) ||
     (maleGenLen.n_elem!=motherGeno.n_elem)){
    Rcpp::stop("genLen must have one value per chromosome");
  }
  if(mother.n_elem>0){
    if((mother.min()<1) || (mother.max()>motherGeno(0).n_slices) ||
       (father.min()<1) || (father.max()>fatherGeno(0).n_slices)){
      Rcpp::stop("Invalid parent index");
    }
  }
  mother -= 1; // R to C++
  father -= 1; // R to C++
  arma::uword ploidy = (motherPloidy+fatherPloidy)/2;
  arma::uword nChr = motherGeno.n_elem;
  arma::uword nInd = mother.n_elem;
  //Output data
  // Sized up front so that gametes are written straight into the output
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  for(arma::uword chr=0; chr<nChr; ++chr){
    geno(chr).set_size(motherGeno(chr).n_rows,ploidy,nInd);
  }
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd,nChr,ploidy);
  }
  arma::uword nBlocks = countBlocks(nInd);
  arma::uword nWork = nChr*nBlocks;
  if(nWork < static_cast<arma::uword>(nThreads) ){
    nThreads = static_cast<int>(nWork);
  }
  if(nThreads < 1){
    nThreads = 1;
  }
  std::vector<alphasimrRng::rngPtr> workRngs = makeWorkRngs(nWork);
  //Loop through chromosome by individual block pairs
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword work=0; work<nWork; ++work){
    arma::uword chr = work/nBlocks;
    arma::uword block = work%nBlocks;
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    alphasimrRng::rngEngine& rng = *workRngs[work];
    arma::Mat<int> hist1, hist2;
    arma::uvec xm(motherPloidy); // Indicator for mother chromosomes
    for(arma::uword i=0; i<motherPloidy; ++i)
      xm(i) = i;
    arma::uvec xf(fatherPloidy); // Indicator for father chromosomes
    for(arma::uword i=0; i<fatherPloidy; ++i)
      xf(i) = i;
    arma::uword progenyChr;
    arma::uword nBins = motherGeno(chr).n_rows;
    arma::Cube<unsigned char>& tmpGeno = geno(chr);

    //Loop through the individuals of this block
    for(arma::uword ind=indStart; ind<indEnd; ++ind){
      progenyChr=0;
      alphasimrRng::shuffle(xm, rng);
      
      //Female gamete
      for(arma::uword x=0; x<motherPloidy; x+=4){
        if((motherPloidy-x)>2){
          if(alphasimrRng::runif(rng)>quadProb){
            //Bivalent 1
            bivalent(motherGeno(chr).slice_colptr(mother(ind), xm(x)),
                     motherGeno(chr).slice_colptr(mother(ind), xm(x+1)),
                     nBins,
                     femaleMap(chr),
                     femaleGenLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x))+1);
              hist1.col(0).replace(200,int(xm(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
            
            //Bivalent 2
            bivalent(motherGeno(chr).slice_colptr(mother(ind), xm(x+2)),
                     motherGeno(chr).slice_colptr(mother(ind), xm(x+3)),
                     nBins,
                     femaleMap(chr),
                     femaleGenLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x+2))+1);
              hist1.col(0).replace(200,int(xm(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
          }else{
            //Quadrivalent
            quadrivalent(motherGeno(chr).slice_colptr(mother(ind), xm(x)),
                         motherGeno(chr).slice_colptr(mother(ind), xm(x+1)),
                         motherGeno(chr).slice_colptr(mother(ind), xm(x+2)),
                         motherGeno(chr).slice_colptr(mother(ind), xm(x+3)),
                         nBins,
                         femaleMap(chr),
                         femaleGenLen(chr),
                         motherCentromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              // Both quadrivalent gametes are labeled 1-4 across all four
              // chromosomes, so each history needs all four mapped back
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x))+1);
              hist1.col(0).replace(200,int(xm(x+1))+1);
              hist1.col(0).replace(300,int(xm(x+2))+1);
              hist1.col(0).replace(400,int(xm(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(xm(x))+1);
              hist2.col(0).replace(200,int(xm(x+1))+1);
              hist2.col(0).replace(300,int(xm(x+2))+1);
              hist2.col(0).replace(400,int(xm(x+3))+1);
              hist.addHist(hist2,ind,chr,progenyChr+1);
            }
            progenyChr += 2;
          }
        }else{
          //Bivalent
          bivalent(motherGeno(chr).slice_colptr(mother(ind), xm(x)),
                   motherGeno(chr).slice_colptr(mother(ind), xm(x+1)),
                   nBins,
                   femaleMap(chr),
                   femaleGenLen(chr),
                   v,
                   p,
                   tmpGeno.slice_colptr(ind, progenyChr),
                   hist1,
                   rng);
          if(trackRec){
            hist1.col(0) *= 100; //To avoid conflicts
            hist1.col(0).replace(100,int(xm(x))+1);
            hist1.col(0).replace(200,int(xm(x+1))+1);
            hist.addHist(hist1,ind,chr,progenyChr);
          }
          ++progenyChr;
        }
      }
      
      //Male gamete
      alphasimrRng::shuffle(xf, rng);
      for(arma::uword x=0; x<fatherPloidy; x+=4){
        if((fatherPloidy-x)>2){
          if(alphasimrRng::runif(rng)>quadProb){
            //Bivalent 1
            bivalent(fatherGeno(chr).slice_colptr(father(ind), xf(x)),
                     fatherGeno(chr).slice_colptr(father(ind), xf(x+1)),
                     nBins,
                     maleMap(chr),
                     maleGenLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x))+1);
              hist1.col(0).replace(200,int(xf(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
            
            //Bivalent 2
            bivalent(fatherGeno(chr).slice_colptr(father(ind), xf(x+2)),
                     fatherGeno(chr).slice_colptr(father(ind), xf(x+3)),
                     nBins,
                     maleMap(chr),
                     maleGenLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x+2))+1);
              hist1.col(0).replace(200,int(xf(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
          }else{
            //Quadrivalent
            quadrivalent(fatherGeno(chr).slice_colptr(father(ind), xf(x)),
                         fatherGeno(chr).slice_colptr(father(ind), xf(x+1)),
                         fatherGeno(chr).slice_colptr(father(ind), xf(x+2)),
                         fatherGeno(chr).slice_colptr(father(ind), xf(x+3)),
                         nBins,
                         maleMap(chr),
                         maleGenLen(chr),
                         fatherCentromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              // Both quadrivalent gametes are labeled 1-4 across all four
              // chromosomes, so each history needs all four mapped back
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x))+1);
              hist1.col(0).replace(200,int(xf(x+1))+1);
              hist1.col(0).replace(300,int(xf(x+2))+1);
              hist1.col(0).replace(400,int(xf(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(xf(x))+1);
              hist2.col(0).replace(200,int(xf(x+1))+1);
              hist2.col(0).replace(300,int(xf(x+2))+1);
              hist2.col(0).replace(400,int(xf(x+3))+1);
              hist.addHist(hist2,ind,chr,progenyChr+1);
            }
            progenyChr += 2;
          }
        }else{
          //Bivalent
          bivalent(fatherGeno(chr).slice_colptr(father(ind), xf(x)),
                   fatherGeno(chr).slice_colptr(father(ind), xf(x+1)),
                   nBins,
                   maleMap(chr),
                   maleGenLen(chr),
                   v,
                   p,
                   tmpGeno.slice_colptr(ind, progenyChr),
                   hist1,
                   rng);
          if(trackRec){
            hist1.col(0) *= 100; //To avoid conflicts
            hist1.col(0).replace(100,int(xf(x))+1);
            hist1.col(0).replace(200,int(xf(x+1))+1);
            hist.addHist(hist1,ind,chr,progenyChr);
          }
          ++progenyChr;
        }
      }
    } //End individual loop
  } //End work loop
  if(trackRec){
    return Rcpp::List::create(Rcpp::Named("geno")=geno,
                              Rcpp::Named("recHist")=hist.pack());
  }
  return Rcpp::List::create(Rcpp::Named("geno")=geno);
}

// Creates DH lines from diploid individuals
// [[Rcpp::export]]
Rcpp::List createDH2(
    const arma::field<arma::Cube<unsigned char> >& geno, 
    const arma::uvec& nDH, const arma::field<arma::vec>& genMap, 
    const arma::vec& genLen, double v, double p, bool trackRec,
    int nThreads){
  arma::uword nChr = geno.n_elem;
  arma::uword nInd = geno(0).n_slices;
  // Each chromosome's genetic length is read inside the parallel region
  if(genLen.n_elem != nChr){
    Rcpp::stop("genLen must have one value per chromosome");
  }
  // nDH holds one value per individual, so an individual's output slices
  // no longer sit at a fixed stride. They start at the running total of
  // the individuals before it. The offsets are built once, before the
  // parallel region, so that every thread can work out its own write
  // positions without coordinating with any other. Accumulating them as
  // the loop ran would make the result depend on the thread count.
  if(nDH.n_elem != nInd){
    Rcpp::stop("nDH must have one value per individual");
  }
  arma::uword nOut = arma::accu(nDH);
  arma::uvec offset = arma::cumsum(nDH) - nDH;
  //Output data
  // Sized up front so that gametes are written straight into the output
  arma::field<arma::Cube<unsigned char> > output(nChr);
  for(arma::uword chr=0; chr<nChr; ++chr){
    output(chr).set_size(geno(chr).n_rows,2,nOut);
  }
  RecHist hist;
  if(trackRec){
    hist.setSize(nOut,nChr,2);
  }
  arma::uword nBlocks = countBlocks(nInd);
  arma::uword nWork = nChr*nBlocks;
  if(nWork < static_cast<arma::uword>(nThreads) ){
    nThreads = static_cast<int>(nWork);
  }
  if(nThreads < 1){
    nThreads = 1;
  }
  std::vector<alphasimrRng::rngPtr> workRngs = makeWorkRngs(nWork);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword work=0; work<nWork; ++work){ //Work loop
    arma::uword chr = work/nBlocks;
    arma::uword block = work%nBlocks;
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    alphasimrRng::rngEngine& rng = *workRngs[work];
    arma::Mat<int> histMat;
    arma::uword nBins = geno(chr).n_rows;
    arma::Cube<unsigned char>& tmp = output(chr);
    arma::uvec x = {0,1};
    for(arma::uword ind=indStart; ind<indEnd; ++ind){ //Individual loop
      for(arma::uword i=0; i<nDH(ind); ++i){ //nDH loop
        arma::uword out = offset(ind)+i;
        alphasimrRng::shuffle(x, rng);
        bivalent(geno(chr).slice_colptr(ind, x(0)),
                 geno(chr).slice_colptr(ind, x(1)),
                 nBins,
                 genMap(chr),
                 genLen(chr),
                 v,
                 p,
                 tmp.slice_colptr(out, 0),
                 histMat,
                 rng);
        // Both haplotypes of a doubled haploid are the same gamete
        std::memcpy(tmp.slice_colptr(out, 1),
                    tmp.slice_colptr(out, 0), nBins);
        for(arma::uword j=0; j<2; ++j){ //ploidy loop
          if(trackRec){
            if((x(0)==1) & (j==0)){
              histMat.col(0).transform([](int val){return val%2+1;});
            }
            hist.addHist(histMat,out,chr,j);
          }
        } //End ploidy loop
      } //End nDH loop
    } //End individual loop
  } //End work loop
  if(trackRec){
    return Rcpp::List::create(Rcpp::Named("geno")=output,
                              Rcpp::Named("recHist")=hist.pack());
  }
  return Rcpp::List::create(Rcpp::Named("geno")=output);
}

// Samples gametes from individuals
// [[Rcpp::export]]
Rcpp::List createReducedGenome(
    const arma::field<arma::Cube<unsigned char> >& geno, 
    const arma::uvec& nProgeny, const arma::field<arma::vec>& genMap, 
    const arma::vec& genLen, double v, double p, bool trackRec,
    arma::uword ploidy, arma::vec& centromere, double quadProb,
    int nThreads){
  arma::uword nChr = geno.n_elem;
  arma::uword nInd = geno(0).n_slices;
  // Each chromosome's genetic length is read inside the parallel region
  if(genLen.n_elem != nChr){
    Rcpp::stop("genLen must have one value per chromosome");
  }
  // nProgeny holds one value per individual. Unlike createDH2 the work
  // loop below is partitioned over progeny rather than over individuals,
  // so each output slice has to be able to name its own parent. With a
  // single value that was ind/nProgeny; with a vector it is a lookup
  // table, built once before the parallel region so that the answer does
  // not depend on how the work was divided.
  if(nProgeny.n_elem != nInd){
    Rcpp::stop("nProgeny must have one value per individual");
  }
  arma::uword nOut = arma::accu(nProgeny);
  arma::uvec parentOf(nOut);
  {
    arma::uword k=0;
    for(arma::uword i=0; i<nInd; ++i){
      for(arma::uword j=0; j<nProgeny(i); ++j){
        parentOf(k) = i;
        ++k;
      }
    }
  }
  //Output data
  // Sized up front so that gametes are written straight into the output
  arma::field<arma::Cube<unsigned char> > output(nChr);
  for(arma::uword chr=0; chr<nChr; ++chr){
    output(chr).set_size(geno(chr).n_rows,ploidy/2,nOut);
  }
  RecHist hist;
  if(trackRec){
    hist.setSize(nOut,nChr,ploidy/2);
  }
  arma::uword nBlocks = countBlocks(nOut);
  arma::uword nWork = nChr*nBlocks;
  if(nWork < static_cast<arma::uword>(nThreads) ){
    nThreads = static_cast<int>(nWork);
  }
  if(nThreads < 1){
    nThreads = 1;
  }
  std::vector<alphasimrRng::rngPtr> workRngs = makeWorkRngs(nWork);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword work=0; work<nWork; ++work){ //Work loop
    arma::uword chr = work/nBlocks;
    arma::uword block = work%nBlocks;
    arma::uword indStart = blockStart(nOut, nBlocks, block);
    arma::uword indEnd = blockStart(nOut, nBlocks, block+1);
    alphasimrRng::rngEngine& rng = *workRngs[work];
    arma::Mat<int> hist1, hist2;
    arma::uword nBins = geno(chr).n_rows;
    arma::Cube<unsigned char>& tmpGeno = output(chr);
    arma::uvec x(ploidy);
    for(arma::uword i=0; i<ploidy; ++i) 
      x(i) = i;
    for(arma::uword ind=indStart; ind<indEnd; ++ind){ //Individual loop
      alphasimrRng::shuffle(x, rng);
      arma::uword progenyChr=0;
      arma::uword par = parentOf(ind);
      for(arma::uword y=0; y<ploidy; y+=4){
        if((ploidy-y)>2){
          if(alphasimrRng::runif(rng)>quadProb){
            //Bivalent 1
            bivalent(geno(chr).slice_colptr(par, x(y)),
                     geno(chr).slice_colptr(par, x(y+1)),
                     nBins,
                     genMap(chr),
                     genLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(x(y))+1);
              hist1.col(0).replace(200,int(x(y+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
            
            //Bivalent 2
            bivalent(geno(chr).slice_colptr(par, x(y+2)),
                     geno(chr).slice_colptr(par, x(y+3)),
                     nBins,
                     genMap(chr),
                     genLen(chr),
                     v,
                     p,
                     tmpGeno.slice_colptr(ind, progenyChr),
                     hist1,
                     rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(x(y+2))+1);
              hist1.col(0).replace(200,int(x(y+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
          }else{
            //Quadrivalent
            quadrivalent(geno(chr).slice_colptr(par, x(y)),
                         geno(chr).slice_colptr(par, x(y+1)),
                         geno(chr).slice_colptr(par, x(y+2)),
                         geno(chr).slice_colptr(par, x(y+3)),
                         nBins,
                         genMap(chr),
                         genLen(chr),
                         centromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              // Both quadrivalent gametes are labeled 1-4 across all four
              // chromosomes, so each history needs all four mapped back
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(x(y))+1);
              hist1.col(0).replace(200,int(x(y+1))+1);
              hist1.col(0).replace(300,int(x(y+2))+1);
              hist1.col(0).replace(400,int(x(y+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(x(y))+1);
              hist2.col(0).replace(200,int(x(y+1))+1);
              hist2.col(0).replace(300,int(x(y+2))+1);
              hist2.col(0).replace(400,int(x(y+3))+1);
              hist.addHist(hist2,ind,chr,progenyChr+1);
            }
            progenyChr += 2;
          }
        }else{
          //Bivalent
          bivalent(geno(chr).slice_colptr(par, x(y)),
                   geno(chr).slice_colptr(par, x(y+1)),
                   nBins,
                   genMap(chr),
                   genLen(chr),
                   v,
                   p,
                   tmpGeno.slice_colptr(ind, progenyChr),
                   hist1,
                   rng);
          if(trackRec){
            hist1.col(0) *= 100; //To avoid conflicts
            hist1.col(0).replace(100,int(x(y))+1);
            hist1.col(0).replace(200,int(x(y+1))+1);
            hist.addHist(hist1,ind,chr,progenyChr);
          }
          ++progenyChr;
        }
      } // End ploidy loop
    } // End individual loop
  } //End work loop
  if(trackRec){
    return Rcpp::List::create(Rcpp::Named("geno")=output,
                              Rcpp::Named("recHist")=hist.pack());
  }
  return Rcpp::List::create(Rcpp::Named("geno")=output);
}
