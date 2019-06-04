// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

class RecHist{
public:
  arma::field< //individual
    arma::field< //chromosome
      arma::field< //ploidy
        arma::Mat<int> > > > hist; //(chr, site)
  void setSize(arma::uword nInd, 
               arma::uword nChr, 
               arma::uword ploidy);
  void addHist(arma::Mat<int>& input, 
               arma::uword nInd, 
               arma::uword chrGroup,
               arma::uword chrInd);
  arma::Mat<int> getHist(arma::uword ind, 
                         arma::uword chr,
                         arma::uword par);
};

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

void RecHist::addHist(arma::Mat<int>& input, 
             arma::uword nInd, 
             arma::uword chrGroup,
             arma::uword chrInd){
  hist(nInd)(chrGroup)(chrInd) = input;
}

arma::Mat<int> RecHist::getHist(arma::uword ind, 
                                arma::uword chr,
                                arma::uword par){
  return hist(ind)(chr)(par);
}

// Samples position of chiasmata following a gamma model
arma::vec sampleChiasmata(double end, double v=1,
                          double start=-10){
  arma::uword n = 35;
  arma::vec output = arma::randg<arma::vec>(n, arma::distr_param(v/2,1/v));
  output = cumsum(output)+start;
  // Add additional values if less than end
  while(output(output.n_elem-1)<end){
    arma::vec tmp = arma::randg<arma::vec>(n, arma::distr_param(v/2,1/v));
    tmp = cumsum(tmp) + output(output.n_elem-1);
    output = join_cols(output, tmp);
  }
  // Remove values less than 0
  output = output(find(output>0));
  // Return values less than the end
  return output(find(output<end));
}

// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns -1 if value is smaller than the values of x
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
  // Works backwards, because the last crossover is observed
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
arma::Mat<int> findBivalentCO(const arma::vec& genMap, 
                              double v){
  arma::uword startPos=0, endPos, readChr=0, nCO;
  double genLen = genMap(genMap.n_elem-1);
  
  // Find crossover positions
  arma::vec posCO = sampleChiasmata(genLen, v);
  if(posCO.n_elem==0){
   arma::Mat<int> output(1,2,arma::fill::ones);
   return output;
  }
  
  // Thin crossovers
  arma::vec thin(posCO.n_elem, arma::fill::randu);
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

// Finds recombination map for cross type quadrivalent
// Chromosome pairing in the head: 1:2 and 3:4
// Chromosome pairing in the tail: 1:4 and 2:3
// chr is the selected chromosome (1-4)
// exchange is the exchange point between chromosomes 
arma::Mat<int>  findQuadrivalentCO(arma::uword chr, //1-4
                                   double exchange, //0-genLen
                                   double centromere, //0-genLen
                                   const arma::vec& genMap, 
                                   double v){ //Ordered values 0-genLen, length nSites
  arma::uvec relPosCO(3,arma::fill::zeros), pairChr1(2), pairChr2(2);
  arma::uword startPos=0, endPos, readChr=0, nCO_1, nCO;
  double genLen = genMap(genMap.n_elem-1);
  
  // Find crossover positions
  arma::vec posCO = sampleChiasmata(genLen, v);
  if(posCO.n_elem==0){
    arma::Mat<int> output(1,2,arma::fill::ones);
    return output;
  }
  
  // Thin crossovers
  arma::vec thin(posCO.n_elem, arma::fill::randu);
  posCO = posCO(find(thin>0.5));
  nCO = posCO.n_elem;
  
  arma::Mat<int> output(nCO+1,2);
  if(nCO==0){
    output(0,0) = chr;
    output(0,1) = 1;
    return output;
  }
  
  //Find starting chromosome pair and number of CO prior to exchange
  //Find second chromosome pair and number of CO after exchange
  if(centromere<exchange){ //Centromere in head
    for(arma::uword i=0; i<nCO; ++i){
      if(posCO(i)<centromere){
        ++relPosCO(0);
      }else if(posCO(i)<exchange){
        ++relPosCO(1);
      }else{
        ++relPosCO(2);
      }
    }
    nCO_1 = relPosCO(0) + relPosCO(1);
    switch(chr){ // Identity of selected centromere
    case 1: 
      if(relPosCO(0)%2 == 0){
        pairChr1(0) = 1;
        pairChr1(1) = 2;
      }else{
        pairChr1(0) = 2;
        pairChr1(1) = 1;
      }
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 1;
        pairChr2(1) = 4;
      }else{
        pairChr2(0) = 2;
        pairChr2(1) = 3;
      }
      break; 
    case 2:
      if(relPosCO(0)%2 == 0){
        pairChr1(0) = 2;
        pairChr1(1) = 1;
      }else{
        pairChr1(0) = 1;
        pairChr1(1) = 2;
      }
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 2;
        pairChr2(1) = 3;
      }else{
        pairChr2(0) = 1;
        pairChr2(1) = 4;
      }
      break; 
    case 3:
      if(relPosCO(0)%2 == 0){
        pairChr1(0) = 3;
        pairChr1(1) = 4;
      }else{
        pairChr1(0) = 4;
        pairChr1(1) = 3;
      }
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 3;
        pairChr2(1) = 2;
      }else{
        pairChr2(0) = 4;
        pairChr2(1) = 1;
      }
      break; 
    case 4:
      if(relPosCO(0)%2 == 0){
        pairChr1(0) = 4;
        pairChr1(1) = 3;
      }else{
        pairChr1(0) = 3;
        pairChr1(1) = 4;
      }
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 4;
        pairChr2(1) = 1;
      }else{
        pairChr2(0) = 3;
        pairChr2(1) = 2;
      }
    }
  }else{ //Centromere in tail
    for(arma::uword i=0; i<nCO; ++i){
      if(posCO(i)<exchange){
        ++relPosCO(0);
      }else if(posCO(i)<centromere){
        ++relPosCO(1);
      }else{
        ++relPosCO(2);
      }
    }
    nCO_1 = relPosCO(0);
    switch(chr){ // Identity of selected centromere
    case 1:
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 1;
        pairChr2(1) = 4;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 1;
          pairChr1(1) = 2;
        }else{
          pairChr1(0) = 2;
          pairChr1(1) = 1;
        }
      }else{
        pairChr2(0) = 4;
        pairChr2(1) = 1;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 4;
          pairChr1(1) = 3;
        }else{
          pairChr1(0) = 3;
          pairChr1(1) = 4;
        }
      }
      break; 
    case 2:
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 2;
        pairChr2(1) = 3;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 2;
          pairChr1(1) = 1;
        }else{
          pairChr1(0) = 1;
          pairChr1(1) = 2;
        }
      }else{
        pairChr2(0) = 3;
        pairChr2(1) = 2;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 3;
          pairChr1(1) = 4;
        }else{
          pairChr1(0) = 4;
          pairChr1(1) = 3;
        }
      }
      break; 
    case 3:
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 3;
        pairChr2(1) = 2;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 3;
          pairChr1(1) = 4;
        }else{
          pairChr1(0) = 4;
          pairChr1(1) = 3;
        }
      }else{
        pairChr2(0) = 2;
        pairChr2(1) = 3;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 2;
          pairChr1(1) = 1;
        }else{
          pairChr1(0) = 1;
          pairChr1(1) = 2;
        }
      }
      break; 
    case 4:
      if(relPosCO(1)%2 == 0){
        pairChr2(0) = 4;
        pairChr2(1) = 1;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 4;
          pairChr1(1) = 3;
        }else{
          pairChr1(0) = 3;
          pairChr1(1) = 4;
        }
      }else{
        pairChr2(0) = 1;
        pairChr2(1) = 4;
        if(relPosCO(0)%2 == 0){
          pairChr1(0) = 1;
          pairChr1(1) = 2;
        }else{
          pairChr1(0) = 2;
          pairChr1(1) = 1;
        }
      }
    }
  }
  
  // Find crossover sites on map
  output(0,0) = pairChr1(0);
  output(0,1) = 1;
  
  for(arma::uword i=0; i<nCO_1; ++i){
    ++readChr;
    readChr = readChr%2;
    endPos = intervalSearch(genMap,posCO(i),startPos);
    output(i+1,0) = pairChr1(readChr);
    output(i+1,1) = endPos+2;
    startPos = endPos;
  }
  readChr = 0;
  for(arma::uword i=nCO_1; i<nCO; ++i){
    ++readChr;
    readChr = readChr%2;
    endPos = intervalSearch(genMap,posCO(i),startPos);
    output(i+1,0) = pairChr2(readChr);
    output(i+1,1) = endPos+2;
    startPos = endPos;
  }
  return removeDoubleCO(output);
}

void transferGeno(const arma::Col<unsigned char>& inChr,
                  arma::Col<unsigned char>& outChr,
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
    inBits = toBits(inChr(startByte));
    outBits = toBits(outChr(startByte));
    if(stopByte > startByte){
      // Transferring more than this byte
      for(int i=startBit; i<8; ++i){
        outBits[i] = inBits[i];
      }
      outChr(startByte) = toByte(outBits);
      startBit = 0;
      ++startByte;
    }else{
      // Only transferring within this byte
      for(int i=startBit; i<stopBit; ++i){
        outBits[i] = inBits[i];
      }
      outChr(startByte) = toByte(outBits);
      return;
    }
  }
  // Transfer full bytes
  if(stopByte > startByte){
    outChr(arma::span(startByte,stopByte-1)) = 
      inChr(arma::span(startByte,stopByte-1));
    startByte = stopByte;
  }
  // Transfer partial stop
  if(inChr.n_elem == startByte){
    // End has been reached
    return;
  }else{
    if(stopBit > startBit){
      inBits = toBits(inChr(startByte));
      outBits = toBits(outChr(startByte));
      for(int i = startBit; i<stopBit; ++i){
        outBits[i] = inBits[i];
      }
      outChr(startByte) = toByte(outBits);
    }
  }
}

//Simulates a gamete using a count-location model for recombination
void bivalent(const arma::Col<unsigned char>& chr1,
              const arma::Col<unsigned char>& chr2,
              const arma::vec& genMap,
              double v,
              arma::Col<unsigned char>& output,
              arma::Mat<int>& hist){
  hist = findBivalentCO(genMap, v);
  if(hist.n_rows==1){
    output = chr1;
  }else{
    int nBins = chr1.n_elem;
    // Fill-in based on recombination history
    for(arma::uword i=0; i<(hist.n_rows-1); ++i){
      switch(hist(i,0)){
      case 1: //Chromosome 1
        transferGeno(chr1, output, 
                     hist(i,1), hist(i+1,1));
        break; 
      case 2: //Chromosome 2
        transferGeno(chr2, output, 
                     hist(i,1), hist(i+1,1));
      }
    }
    // Fill-in last sites
    switch(hist(hist.n_rows-1,0)){
    case 1:
      transferGeno(chr1, output, 
                   hist(hist.n_rows-1,1), 
                   nBins*8+1);
      break; 
    case 2:
      transferGeno(chr2, output, 
                   hist(hist.n_rows-1,1), 
                   nBins*8+1);
    }
  }
}

//Simulates a gamete using a count-location model for recombination
void quadrivalent(const arma::Col<unsigned char>& chr1,
                  const arma::Col<unsigned char>& chr2,
                  const arma::Col<unsigned char>& chr3,
                  const arma::Col<unsigned char>& chr4,
                  const arma::vec& genMap,
                  double v,
                  double centromere,
                  arma::Col<unsigned char>& output1,
                  arma::Col<unsigned char>& output2,
                  arma::Mat<int>& hist1,
                  arma::Mat<int>& hist2){
  int nBins = chr1.n_elem;
  
  //Find exchange point
  double genLen = genMap(genMap.n_elem-1);
  arma::vec exchange(1,arma::fill::randu);
  exchange *= genLen;
  
  //Sample chromosomes
  arma::uvec selChr = sampleInt(2,4);
  ++selChr; //1-4 coding
  
  //Resolve first gamete
  hist1 = findQuadrivalentCO(selChr(0), exchange(0),
                             centromere, genMap, v);
  if(hist1.n_rows==1){
    switch(hist1(0,0)){
    case 1:
      output1 = chr1;
      break;
    case 2:
      output1 = chr2;
      break;
    case 3:
      output1 = chr3;
      break;
    case 4:
      output1 = chr4;
    }
  }else{
    // Fill-in based on recombination history
    for(arma::uword i=0; i<(hist1.n_rows-1); ++i){
      switch(hist1(i,0)){
      case 1:
        transferGeno(chr1, output1, 
                     hist1(i,1), hist1(i+1,1));
        break;  
      case 2:
        transferGeno(chr2, output1, 
                     hist1(i,1), hist1(i+1,1));
        break;   
      case 3:
        transferGeno(chr3, output1, 
                     hist1(i,1), hist1(i+1,1));
        break;  
      case 4:
        transferGeno(chr4, output1, 
                     hist1(i,1), hist1(i+1,1));
      }
    }
    // Fill-in last sites
    switch(hist1(hist1.n_rows-1,0)){
    case 1:
      transferGeno(chr1, output1, 
                   hist1(hist1.n_rows-1,1), 
                   nBins*8+1);
      break;
    case 2:
      transferGeno(chr2, output1, 
                   hist1(hist1.n_rows-1,1), 
                   nBins*8+1);
      break;
    case 3:
      transferGeno(chr3, output1, 
                   hist1(hist1.n_rows-1,1), 
                   nBins*8+1);
      break;
    case 4:
      transferGeno(chr4, output1, 
                   hist1(hist1.n_rows-1,1), 
                   nBins*8+1);
    }
  }
  
  //Resolve second gamete
  hist2 = findQuadrivalentCO(selChr(1), exchange(0),
                             centromere, genMap, v);
  if(hist2.n_rows==1){
    switch(hist2(0,0)){
    case 1:
      output2 = chr1;
      break;
    case 2:
      output2 = chr2;
      break;
    case 3:
      output2 = chr3;
      break;
    case 4:
      output2 = chr4;
    }
  }else{
    // Fill-in based on recombination history
    for(arma::uword i=0; i<(hist2.n_rows-1); ++i){
      switch(hist2(i,0)){
      case 1:
        transferGeno(chr1, output2, 
                     hist2(i,1), hist2(i+1,1));
        break;  
      case 2:
        transferGeno(chr2, output2, 
                     hist2(i,1), hist2(i+1,1));
        break; 
      case 3:
        transferGeno(chr3, output2, 
                     hist2(i,1), hist2(i+1,1));
        break;  
      case 4:
        transferGeno(chr4, output2, 
                     hist2(i,1), hist2(i+1,1));
      }
    }
    // Fill-in last sites
    switch(hist2(hist2.n_rows-1,0)){
    case 1:
      transferGeno(chr1, output2, 
                   hist2(hist2.n_rows-1,1), 
                   nBins*8+1);
      break; 
    case 2:
      transferGeno(chr2, output2, 
                   hist2(hist2.n_rows-1,1), 
                   nBins*8+1);
      break; 
    case 3:
      transferGeno(chr3, output2, 
                   hist2(hist2.n_rows-1,1), 
                   nBins*8+1);
      break; 
    case 4:
      transferGeno(chr4, output2, 
                   hist2(hist2.n_rows-1,1), 
                   nBins*8+1);
    }
  }
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
// v: interference paramater for gamma model
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
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd,nChr,ploidy);
  }
  //Loop through chromosomes
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword chr=0; chr<nChr; ++chr){
    arma::vec u(1);
    arma::Mat<int> hist1, hist2;
    arma::uvec xm(motherPloidy); // Indicator for mother chromosomes
    for(arma::uword i=0; i<motherPloidy; ++i)
      xm(i) = i;
    arma::uvec xf(fatherPloidy); // Indicator for father chromosomes
    for(arma::uword i=0; i<fatherPloidy; ++i)
      xf(i) = i;
    arma::uword progenyChr;
    arma::uword nBins = motherGeno(chr).n_rows;
    arma::Cube<unsigned char> tmpGeno(nBins,ploidy,nInd);
    arma::Col<unsigned char> gamete1(nBins), gamete2(nBins);
    //Loop through individuals
    for(arma::uword ind=0; ind<nInd; ++ind){
      progenyChr=0;
      xm = shuffle(xm);
      //Female gamete
      for(arma::uword x=0; x<motherPloidy; x+=4){
        if((motherPloidy-x)>2){
          u.randu();
          if(u(0)>quadProb){
            //Bivalent 1
            bivalent(motherGeno(chr).slice(mother(ind)).col(xm(x)),
                     motherGeno(chr).slice(mother(ind)).col(xm(x+1)),
                     femaleMap(chr),
                     v,
                     gamete1,
                     hist1);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x))+1);
              hist1.col(0).replace(200,int(xm(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
            //Bivalent 2
            bivalent(motherGeno(chr).slice(mother(ind)).col(xm(x+2)),
                     motherGeno(chr).slice(mother(ind)).col(xm(x+3)),
                     femaleMap(chr),
                     v,
                     gamete1,
                     hist1);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xm(x+2))+1);
              hist1.col(0).replace(200,int(xm(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
          }else{
            //Quadrivalent
            quadrivalent(motherGeno(chr).slice(mother(ind)).col(xm(x)),
                         motherGeno(chr).slice(mother(ind)).col(xm(x+1)),
                         motherGeno(chr).slice(mother(ind)).col(xm(x+2)),
                         motherGeno(chr).slice(mother(ind)).col(xm(x+3)),
                         femaleMap(chr),
                         v,
                         motherCentromere(chr),
                         gamete1,
                         gamete2,
                         hist1,
                         hist2);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            tmpGeno.slice(ind).col(progenyChr+1) = gamete2;
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
          bivalent(motherGeno(chr).slice(mother(ind)).col(xm(x)),
                   motherGeno(chr).slice(mother(ind)).col(xm(x+1)),
                   femaleMap(chr),
                   v,
                   gamete1,
                   hist1);
          tmpGeno.slice(ind).col(progenyChr) = gamete1;
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
      xf = shuffle(xf);
      for(arma::uword x=0; x<fatherPloidy; x+=4){
        if((fatherPloidy-x)>2){
          u.randu();
          if(u(0)>quadProb){
            //Bivalent 1
            bivalent(fatherGeno(chr).slice(father(ind)).col(xf(x)),
                     fatherGeno(chr).slice(father(ind)).col(xf(x+1)),
                     maleMap(chr),
                     v,
                     gamete1,
                     hist1);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x))+1);
              hist1.col(0).replace(200,int(xf(x+1))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
            //Bivalent 2
            bivalent(fatherGeno(chr).slice(father(ind)).col(xf(x+2)),
                     fatherGeno(chr).slice(father(ind)).col(xf(x+3)),
                     maleMap(chr),
                     v,
                     gamete1,
                     hist1);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            if(trackRec){
              hist1.col(0) *= 100; //To avoid conflicts
              hist1.col(0).replace(100,int(xf(x+2))+1);
              hist1.col(0).replace(200,int(xf(x+3))+1);
              hist.addHist(hist1,ind,chr,progenyChr);
            }
            ++progenyChr;
          }else{
            //Quadrivalent
            quadrivalent(fatherGeno(chr).slice(father(ind)).col(xf(x)),
                         fatherGeno(chr).slice(father(ind)).col(xf(x+1)),
                         fatherGeno(chr).slice(father(ind)).col(xf(x+2)),
                         fatherGeno(chr).slice(father(ind)).col(xf(x+3)),
                         maleMap(chr),
                         v,
                         fatherCentromere(chr),
                         gamete1,
                         gamete2,
                         hist1,
                         hist2);
            tmpGeno.slice(ind).col(progenyChr) = gamete1;
            tmpGeno.slice(ind).col(progenyChr+1) = gamete2;
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
          bivalent(fatherGeno(chr).slice(father(ind)).col(xf(x)),
                   fatherGeno(chr).slice(father(ind)).col(xf(x+1)),
                   maleMap(chr),
                   v,
                   gamete1,
                   hist1);
          tmpGeno.slice(ind).col(progenyChr) = gamete1;
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
    geno(chr) = tmpGeno;
  } //End chromosome loop
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
    double v, bool trackRec, int nThreads){
  arma::uword nChr = geno.n_elem;
  arma::uword nInd = geno(0).n_slices;
  //Output data
  arma::field<arma::Cube<unsigned char> > output(nChr);
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd*nDH,nChr,2);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword chr=0; chr<nChr; ++chr){ //Chromosome loop
    arma::Mat<int> histMat;
    arma::uword nBins = geno(chr).n_rows;
    arma::Cube<unsigned char> tmp(nBins,2,nInd*nDH);
    arma::Col<unsigned char> gamete(nBins);
    arma::uvec x = {0,1};
    for(arma::uword ind=0; ind<nInd; ++ind){ //Individual loop
      for(arma::uword i=0; i<nDH; ++i){ //nDH loop
        x = shuffle(x);
        bivalent(geno(chr).slice(ind).col(x(0)),
                 geno(chr).slice(ind).col(x(1)),
                 genMap(chr),
                 v,
                 gamete,
                 histMat);
        for(arma::uword j=0; j<2; ++j){ //ploidy loop
          tmp.slice(i+ind*nDH).col(j) = gamete;
          if(trackRec){
            if(x(0)==1){
              histMat.col(0).transform([](int val){return val%2+1;});
            }
            hist.addHist(histMat,i+ind*nDH,chr,j);
          }
        } //End ploidy loop
      } //End nDH loop
    } //End individual loop
    output(chr) = tmp;
  } //End chromosome loop
  if(trackRec){
    return Rcpp::List::create(Rcpp::Named("geno")=output,
                              Rcpp::Named("recHist")=hist.hist);
  }
  return Rcpp::List::create(Rcpp::Named("geno")=output);
}

// Converts    recHist (recombinations between generations) to
//          ibdRecHist (recombinations since the base generation/population)
// [[Rcpp::export]]
Rcpp::List getIbdRecHist(const Rcpp::List          & recHist,
                         const Rcpp::IntegerMatrix & pedigree,
                         const Rcpp::IntegerVector & nLociPerChr) {
  // This is an utterly complicated function! There has to be a neater way to do this. Gregor
  RecHist ibdRecHist;
  arma::uword nInd = pedigree.nrow();
  arma::uword nChr = nLociPerChr.size();
  ibdRecHist.setSize(nInd, nChr, 2);
  for (arma::uword ind = 0; ind < nInd; ++ind) {
    Rcpp::List recHistInd = recHist(ind);
    for (arma::uword par = 0; par < 2; ++par) {
      int pId = pedigree(ind, par);
      if (pId == 0) { // Individual is     a founder --> set founder gamete code
        for (arma::uword chr = 0; chr < nChr; ++chr) {
          if (0 < nLociPerChr(chr)) {
            arma::Mat<int> recHistIndChrPar;
            recHistIndChrPar.set_size(1, 2);
            recHistIndChrPar(0, 0) = 2 * (ind + 1) - 1 + par;
            recHistIndChrPar(0, 1) = 1;
            ibdRecHist.addHist(recHistIndChrPar, ind, chr, par);
          }
        }
      } else {        // Individual is not a founder --> get founder gamete code & recombinations
        pId -= 1; // R to C++ indexing
        Rcpp::List recHistPar = recHist(pId);
        for (arma::uword chr = 0; chr < nChr; ++chr) {
          if (0 < nLociPerChr(chr)) {
            Rcpp::List recHistIndChr = recHistInd(chr);
            arma::Mat<int> recHistIndChrPar = recHistIndChr(par);
            arma::uword nRecSegInd = recHistIndChrPar.n_rows;
            if (recHistPar.size() == 0) { // Parent is     a founder and has no recHist info --> get founder gamete codes and put them onto individual recombinations
              for (arma::uword recSegInd = 0; recSegInd < nRecSegInd; ++recSegInd) {
                int source = recHistIndChrPar(recSegInd, 0) - 1;
                recHistIndChrPar(recSegInd, 0) = ibdRecHist.getHist(pId, chr, source)(0, 0);
              }
              ibdRecHist.addHist(recHistIndChrPar, ind, chr, par);
            } else {                      // Parent is not a founder and has    recHist info --> parse and combine parent and individual recombinations
              // Parent's all ancestral recombinations
              arma::Mat<int> ibdRecHistParChrPar1 = ibdRecHist.getHist(pId, chr, 0);
              arma::Mat<int> ibdRecHistParChrPar2 = ibdRecHist.getHist(pId, chr, 1);
              arma::field<arma::Mat<int> > ibdRecHistParChrPar(2);
              ibdRecHistParChrPar(0) = ibdRecHistParChrPar1;
              ibdRecHistParChrPar(1) = ibdRecHistParChrPar2;
              arma::uvec nIbdRecSegParChrPar(2);
              nIbdRecSegParChrPar(0) = ibdRecHistParChrPar(0).n_rows;
              nIbdRecSegParChrPar(1) = ibdRecHistParChrPar(1).n_rows;
              
              // Find and advance the ancestral recombinations in line with the recent (parent-progeny) recombinations
              arma::uvec ibdRecSegPar(2);
              int nIbdSegInd;
              arma::Mat<int> ibdRecHistIndChrPar;
              for (arma::uword run = 0; run < 2; ++run) {
                if (run != 0) {
                  ibdRecHistIndChrPar.set_size(nIbdSegInd, 2);
                }
                ibdRecSegPar(0) = 0;
                ibdRecSegPar(1) = 0;
                nIbdSegInd = 0;
                for (arma::uword recSegInd = 0; recSegInd < nRecSegInd; ++recSegInd) {
                  int source = recHistIndChrPar(recSegInd, 0) - 1;
                  int startInd = recHistIndChrPar(recSegInd, 1);
                  int stopInd;
                  if (recSegInd == (nRecSegInd - 1)) {
                    stopInd = nLociPerChr(chr);
                  } else {
                    stopInd = recHistIndChrPar(recSegInd + 1, 1) - 1;
                  }

                  bool loop = true;
                  while (loop & (ibdRecSegPar(source) < nIbdRecSegParChrPar(source))) {
                    int sourcePar = ibdRecHistParChrPar(source)(ibdRecSegPar(source), 0);
                    int startPar  = ibdRecHistParChrPar(source)(ibdRecSegPar(source), 1);
                    int stopPar;
                    if (ibdRecSegPar(source) == (nIbdRecSegParChrPar(source) - 1)) {
                      stopPar = nLociPerChr(chr);
                    } else {
                      stopPar = ibdRecHistParChrPar(source)(ibdRecSegPar(source) + 1, 1) - 1;
                    }

                    if (startInd <= stopPar) {
                      if (stopInd >= startPar) {
                        int startIbd = std::max(startInd, startPar);
                        if (run == 1) {
                          ibdRecHistIndChrPar(nIbdSegInd, 0) = sourcePar;
                          ibdRecHistIndChrPar(nIbdSegInd, 1) = startIbd;
                        }
                        nIbdSegInd += 1;

                        if (stopInd <= stopPar) {
                          loop = false;
                        }
                        if ((stopInd >= stopPar) & (stopPar < nLociPerChr(chr))) {
                          ibdRecSegPar(source) += 1;
                        }
                      } else {
                        loop = false;
                      }
                    } else {
                      ibdRecSegPar(source) += 1;
                    }
                  }
                }
              }
              ibdRecHist.addHist(ibdRecHistIndChrPar, ind, chr, par);
            }
          }
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("ibdRecHist") = ibdRecHist.hist);
}
