// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

class RecHist{
public:
  arma::field< //individual
    arma::field< //chromosome
      arma::field< //ploidy
        arma::Mat<int> > > > hist; //(chr, site)
  void setSize(arma::uword nInd, arma::uword nChr, 
               arma::uword ploidy);
  void addHist(arma::Mat<int>& input, 
               arma::uword nInd, 
               arma::uword chrGroup,
               arma::uword chrInd);
};
void RecHist::setSize(arma::uword nInd, arma::uword nChr, 
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
  // Eliminate rows with no information
  arma::uvec take = find(input.col(0)>0);
  // Eliminate rows where chomosome doesn't change
  arma::uvec takeTake(take.n_elem,arma::fill::zeros);
  takeTake(0) = 1; // First row has starting chromosome
  int lastChr = input(0,0);
  for(arma::uword i=1; i<take.n_elem; ++i){
   if(lastChr!=input(take(i),0)){
     lastChr = input(take(i),0);
     takeTake(i) = 1;
   }
  }
  take = take(find(takeTake>0));
  hist(nInd)(chrGroup)(chrInd) = input.rows(take);
}


// Searches for an interval in x containing value
// Result reported as left most element of the interval
// Returns -1 if value is smaller than the values of x
// Returns last element if value is greater than values of x
// Set left to the smallest value of the interval to search
int intervalSearch(arma::vec x, double value, int left=0){
  // Check if crossover is before beginning
  if(x[left]>value){
    // Return error
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

//Simulates a gamete using a count-location model for recombination, ploidy=2
arma::Col<unsigned char> bivalent(const arma::Col<unsigned char>& chr1,
                                  const arma::Col<unsigned char>& chr2,
                                  const arma::vec& genMap, 
                                  arma::Mat<int>& hist, bool trackRec){
  int nSites = chr1.n_elem;
  double genLen = genMap(nSites-1);
  arma::Col<unsigned char> gamete(nSites);
  // Sample number of chromosomes
  int nCO = Rcpp::rpois(1, genLen)(0);
  // Randomly pick starting chromosome
  int readChr = Rcpp::rbinom(1,1,0.5)(0);
  // Track starting chromosome
  if(trackRec){
    hist.set_size(nCO+2,2);
    hist.fill(0);
    hist(0,0) = readChr+1;
    hist(0,1) = 1;
  }
  if(nCO==0){
    // No CO
    if(readChr){
      gamete = chr2;
    }else{
      gamete = chr1;
    }
  }else{
    // COs present, Create CO locations
    arma::vec posCO(nCO,arma::fill::randu);
    posCO = posCO*genLen;
    posCO = arma::sort(posCO);
    int startPos = 0;
    int endPos;
    if(readChr){
      gamete(0) = chr2(0);
    }else{
      gamete(0) = chr1(0);
    }
    for(arma::uword i=0; i<nCO; ++i){
      endPos = intervalSearch(genMap,posCO[i],startPos);
      // Fill gamete
      if(endPos>startPos){ // Check for double crossovers
        // Fill in segSites
        if(readChr){
          gamete(arma::span(startPos+1,endPos)) = chr2(arma::span(startPos+1,endPos));
        }else{
          gamete(arma::span(startPos+1,endPos)) = chr1(arma::span(startPos+1,endPos));
        }
        if(trackRec){
          hist(i+1,0) = readChr+1;
          hist(i+1,1) = startPos+2;
        }
      }
      startPos = endPos;
      // Switch chromosome
      readChr = ++readChr % 2;
    }
    // Fill in last segSites if needed
    if(endPos<(nSites-1)){
      if(readChr){
        gamete(arma::span(endPos+1,nSites-1)) = chr2(arma::span(endPos+1,nSites-1));
      }else{
        gamete(arma::span(endPos+1,nSites-1)) = chr1(arma::span(endPos+1,nSites-1));
      }
      if(trackRec){
        hist(nCO+1,0) = readChr+1;
        hist(nCO+1,1) = endPos+2;
      }
    }
  }
  return gamete;
}

// Makes crosses between diploid individuals.
// motherGeno: female genotypes
// mother: female parents
// fatherGeno: male genotypes
// father: male parents
// genMaps: chromosome genetic maps
// [[Rcpp::export]]
Rcpp::List cross2(
    const arma::field<arma::Cube<unsigned char> >& motherGeno, 
    arma::uvec mother,
    const arma::field<arma::Cube<unsigned char> >& fatherGeno, 
    arma::uvec father,
    const arma::field<arma::vec>& genMaps,
    double recombRatio, bool trackRec){
  mother -= 1; // R to C++
  father -= 1; // R to C++
  int nChr = motherGeno.n_elem;
  int nInd = mother.n_elem;
  //Output data
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  double femaleRecRate = 2/(1/recombRatio+1);
  double maleRecRate = 2/(recombRatio+1);
  arma::Mat<int> histMat;
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd,nChr,2);
  }
  //Loop through chromosomes
  for(arma::uword chr=0; chr<nChr; ++chr){
    arma::vec maleMap = maleRecRate*genMaps(chr);
    arma::vec femaleMap = femaleRecRate*genMaps(chr);
    int segSites = motherGeno(chr).n_rows;
    arma::Cube<unsigned char> tmpGeno(segSites,2,nInd);
    //Loop through individuals
    for(arma::uword ind=0; ind<nInd; ++ind){
      //Female gamete
      tmpGeno.slice(ind).col(0) = 
        bivalent(motherGeno(chr).slice(mother(ind)).col(0),
                 motherGeno(chr).slice(mother(ind)).col(1),
                 femaleMap,histMat,trackRec);
      if(trackRec){
        hist.addHist(histMat,ind,chr,0);
      }
      //Male gamete
      tmpGeno.slice(ind).col(1) = 
        bivalent(fatherGeno(chr).slice(father(ind)).col(0),
                 fatherGeno(chr).slice(father(ind)).col(1),
                 maleMap,histMat,trackRec);
      if(trackRec){
        hist.addHist(histMat,ind,chr,1);
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
    int nDH, const arma::field<arma::vec>& genMaps,
    double recombRatio, bool useFemale, bool trackRec){
  int nChr = geno.n_elem;
  int nInd = geno(0).n_slices;
  double ratio;
  if(useFemale){
    ratio = 2/(1/recombRatio+1);
  }else{
    ratio = 2/(recombRatio+1);
  }
  //Output data
  arma::field<arma::Cube<unsigned char> > output(nChr);
  arma::Mat<int> histMat;
  RecHist hist;
  if(trackRec){
    hist.setSize(nInd,nChr,2);
  }
  for(arma::uword chr=0; chr<nChr; ++chr){ //Chromosome loop
    arma::vec genMap = ratio*genMaps(chr);
    int segSites = geno(chr).n_rows;
    arma::Cube<unsigned char> tmp(segSites,2,nInd*nDH);
    for(arma::uword ind=0; ind<nInd; ++ind){ //Individual loop
      for(arma::uword i=0; i<nDH; ++i){ //nDH loop
        arma::Col<unsigned char> gamete = 
          bivalent(geno(chr).slice(ind).col(0),
                   geno(chr).slice(ind).col(1),
                   genMap,histMat,trackRec);
        for(arma::uword j=0; j<2; ++j){ //ploidy loop
          tmp.slice(i+ind*nDH).col(j) = gamete;
          if(trackRec){
            hist.addHist(histMat,ind,chr,j);
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
