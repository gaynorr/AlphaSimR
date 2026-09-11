#include "alphasimr.h"
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
  
  // Access history
  arma::Mat<int> getHist(arma::uword ind, 
                         arma::uword chr,
                         arma::uword par);
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
    arma::vec x = alphasimrRng::runifVec(n, rng);
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


// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns an error if value is smaller than the values of x
// Returns last element if value is greater than values of x
// Set left to the smallest value of the interval to search
arma::uword intervalSearch(const arma::vec& x, double& value, arma::uword left=0){
  // Check if crossover is before beginning
  if(x[left]>value){
    Rcpp::stop("intervalSearch searching in impossible interval");
  }
  arma::uword end = x.n_elem-1;
  
  // Check if crossover is at or past end
  if(x[end]<=value){
    return end;
  }
  
  // Perform search
  arma::uword right = end;
  while((right-left)>1){ // Interval can be decreased
    arma::uword middle = (left + right) / 2;
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

// Removes hidden crossovers from recombination map
// Assumes first row is always site 1 and no other row
// will have a value of 1. This logic is based on 
// the implementation of intervalSearch.
arma::Mat<int> removeDoubleCO(const arma::Mat<int>& X){
  if(X.n_rows<3){
    return X;
  }
  
  // Initially assume all rows are useful
  arma::Col<int> take(X.n_rows,arma::fill::ones);
  
  // Remove unobserved crossovers (site doesn't change)
  // Works backward, because the last crossover is observed
  for(arma::uword i=(X.n_rows-2); i>0; --i){
    if(X(i,1) == X(i+1,1)){
      take(i) = 0;
    }
  }
  
  // Remove redundant records (chromosome doesn't change)
  int lastChr = X(0,0);
  for(arma::uword i=1; i<X.n_rows; ++i){
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
// v, the interference parameter
// p, the proportion of non-interfering crossovers
// rng, the explicit dqrng stream used for all random draws in this call
arma::Mat<int> findBivalentCO(const arma::vec& genMap, double v, double p,
                              alphasimrRng::rngEngine& rng){
  arma::uword startPos=0, endPos, readChr=0, nCO;
  double genLen = genMap(genMap.n_elem-1);
  
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
  output.row(0).ones();
  for(arma::uword i=0; i<nCO; ++i){
    ++readChr;
    readChr = readChr%2;
    endPos = intervalSearch(genMap,posCO(i),startPos);
    output(i+1,0) = readChr+1;
    output(i+1,1) = endPos+2;
    startPos = endPos;
  }
  
  return removeDoubleCO(output);
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
                                                double centromere, double v,
                                                double p,
                                                alphasimrRng::rngEngine& rng){
  arma::field<arma::Mat<int> > output(2);
  double genLen = genMap(genMap.n_elem-1);
  
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
    arma::uword startPos=0, endPos, nCO=0;
    output(i)(0,0) = currentChromosome;
    output(i)(0,1) = 1;
    
    if(currentChromosome<3){ // Start in arm 0
      
      // Fill crossovers in the head
      for(arma::uword j=0; j<posCO(0).n_elem; ++j){
        
        // Check if chiasmata involves current chromatid
        switch(currentChromosome){
        case 1: 
          if(chromatidPairs(0)(j,0) == currentChromatid){
            currentChromosome = 2;
            currentChromatid = chromatidPairs(0)(j,1);
            ++nCO;
            endPos = intervalSearch(genMap,posCO(0)(j),startPos);
            output(i)(nCO,0) = currentChromosome;
            output(i)(nCO,1) = endPos+2;
            startPos = endPos;
          }
          break;
        case 2:
          if(chromatidPairs(0)(j,1) == currentChromatid){
            currentChromosome = 1;
            currentChromatid = chromatidPairs(0)(j,0);
            ++nCO;
            endPos = intervalSearch(genMap,posCO(0)(j),startPos);
            output(i)(nCO,0) = currentChromosome;
            output(i)(nCO,1) = endPos+2;
            startPos = endPos;
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
              ++nCO;
              endPos = intervalSearch(genMap,posCO(3)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
            }
            break;
          case 4:
            if(chromatidPairs(3)(j,1) == currentChromatid){
              currentChromosome = 1;
              currentChromatid = chromatidPairs(3)(j,0);
              ++nCO;
              endPos = intervalSearch(genMap,posCO(3)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
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
              ++nCO;
              endPos = intervalSearch(genMap,posCO(1)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
            }
            break;
          case 3:
            if(chromatidPairs(1)(j,1) == currentChromatid){
              currentChromosome = 2;
              currentChromatid = chromatidPairs(1)(j,0);
              ++nCO;
              endPos = intervalSearch(genMap,posCO(1)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
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
            ++nCO;
            endPos = intervalSearch(genMap,posCO(2)(j),startPos);
            output(i)(nCO,0) = currentChromosome;
            output(i)(nCO,1) = endPos+2;
            startPos = endPos;
          }
          break;
        case 4:
          if(chromatidPairs(2)(j,1) == currentChromatid){
            currentChromosome = 3;
            currentChromatid = chromatidPairs(2)(j,0);
            ++nCO;
            endPos = intervalSearch(genMap,posCO(2)(j),startPos);
            output(i)(nCO,0) = currentChromosome;
            output(i)(nCO,1) = endPos+2;
            startPos = endPos;
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
              ++nCO;
              endPos = intervalSearch(genMap,posCO(3)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
            }
            break;
          case 4:
            if(chromatidPairs(3)(j,1) == currentChromatid){
              currentChromosome = 1;
              currentChromatid = chromatidPairs(3)(j,0);
              ++nCO;
              endPos = intervalSearch(genMap,posCO(3)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
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
              ++nCO;
              endPos = intervalSearch(genMap,posCO(1)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
            }
            break;
          case 3:
            if(chromatidPairs(1)(j,1) == currentChromatid){
              currentChromosome = 2;
              currentChromatid = chromatidPairs(1)(j,0);
              ++nCO;
              endPos = intervalSearch(genMap,posCO(1)(j),startPos);
              output(i)(nCO,0) = currentChromosome;
              output(i)(nCO,1) = endPos+2;
              startPos = endPos;
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
              double v,
              double p,
              unsigned char* output,
              arma::Mat<int>& hist,
              alphasimrRng::rngEngine& rng){
  hist = findBivalentCO(genMap, v, p, rng);
  if(hist.n_rows==1){
    // No crossovers, so the gamete is a copy of chromosome 1
    std::memcpy(output, chr1, nBins);
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
  output = findQuadrivalentCO(genMap, centromere, v, p, rng);

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
    bool trackRec,
    arma::uword motherPloidy,
    arma::uword fatherPloidy,
    double v,
    double p,
    const arma::vec& motherCentromere,
    const arma::vec& fatherCentromere,
    double quadProb,
    int nThreads){
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
                         motherCentromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x))+1);
              hist1.col(0).replace(200,int(xm(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(xm(x+2))+1);
              hist2.col(0).replace(200,int(xm(x+3))+1);
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
                         fatherCentromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x))+1);
              hist1.col(0).replace(200,int(xf(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(xf(x+2))+1);
              hist2.col(0).replace(200,int(xf(x+3))+1);
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
                              Rcpp::Named("recHist")=hist.hist);
  }
  return Rcpp::List::create(Rcpp::Named("geno")=geno);
}

// Creates DH lines from diploid individuals
// [[Rcpp::export]]
Rcpp::List createDH2(
    const arma::field<arma::Cube<unsigned char> >& geno, 
    arma::uword nDH, const arma::field<arma::vec>& genMap, 
    double v, double p, bool trackRec, int nThreads){
  arma::uword nChr = geno.n_elem;
  arma::uword nInd = geno(0).n_slices;
  //Output data
  // Sized up front so that gametes are written straight into the output
  arma::field<arma::Cube<unsigned char> > output(nChr);
  for(arma::uword chr=0; chr<nChr; ++chr){
    output(chr).set_size(geno(chr).n_rows,2,nInd*nDH);
  }
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd*nDH,nChr,2);
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
      for(arma::uword i=0; i<nDH; ++i){ //nDH loop
        alphasimrRng::shuffle(x, rng);
        bivalent(geno(chr).slice_colptr(ind, x(0)),
                 geno(chr).slice_colptr(ind, x(1)),
                 nBins,
                 genMap(chr),
                 v,
                 p,
                 tmp.slice_colptr(i+ind*nDH, 0),
                 histMat,
                 rng);
        // Both haplotypes of a doubled haploid are the same gamete
        std::memcpy(tmp.slice_colptr(i+ind*nDH, 1),
                    tmp.slice_colptr(i+ind*nDH, 0), nBins);
        for(arma::uword j=0; j<2; ++j){ //ploidy loop
          if(trackRec){
            if((x(0)==1) & (j==0)){
              histMat.col(0).transform([](int val){return val%2+1;});
            }
            hist.addHist(histMat,i+ind*nDH,chr,j);
          }
        } //End ploidy loop
      } //End nDH loop
    } //End individual loop
  } //End work loop
  if(trackRec){
    return Rcpp::List::create(Rcpp::Named("geno")=output,
                              Rcpp::Named("recHist")=hist.hist);
  }
  return Rcpp::List::create(Rcpp::Named("geno")=output);
}

// Samples gametes from individuals
// [[Rcpp::export]]
Rcpp::List createReducedGenome(
    const arma::field<arma::Cube<unsigned char> >& geno, 
    arma::uword nProgeny, const arma::field<arma::vec>& genMap, 
    double v, double p, bool trackRec, arma::uword ploidy,  
    arma::vec& centromere, double quadProb, int nThreads){
  arma::uword nChr = geno.n_elem;
  arma::uword nInd = geno(0).n_slices;
  //Output data
  // Sized up front so that gametes are written straight into the output
  arma::field<arma::Cube<unsigned char> > output(nChr);
  for(arma::uword chr=0; chr<nChr; ++chr){
    output(chr).set_size(geno(chr).n_rows,ploidy/2,nInd*nProgeny);
  }
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd*nProgeny,nChr,ploidy/2);
  }
  arma::uword nBlocks = countBlocks(nInd*nProgeny);
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
    arma::uword indStart = blockStart(nInd*nProgeny, nBlocks, block);
    arma::uword indEnd = blockStart(nInd*nProgeny, nBlocks, block+1);
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
      arma::uword par = ind/nProgeny;
      for(arma::uword y=0; y<ploidy; y+=4){
        if((ploidy-y)>2){
          if(alphasimrRng::runif(rng)>quadProb){
            //Bivalent 1
            bivalent(geno(chr).slice_colptr(par, x(y)),
                     geno(chr).slice_colptr(par, x(y+1)),
                     nBins,
                     genMap(chr),
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
                         centromere(chr),
                         v,
                         p,
                         tmpGeno.slice_colptr(ind, progenyChr),
                         tmpGeno.slice_colptr(ind, progenyChr+1),
                         hist1,
                         hist2,
                         rng);
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(x(y))+1);
              hist1.col(0).replace(200,int(x(y+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
              hist2.col(0) *= 100; //To avoid conflicts
              hist2.col(0).replace(100,int(x(y+2))+1);
              hist2.col(0).replace(200,int(x(y+3))+1);
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
                              Rcpp::Named("recHist")=hist.hist);
  }
  return Rcpp::List::create(Rcpp::Named("geno")=output);
}
