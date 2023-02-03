#include "alphasimr.h"

// New class for tuning an AD trait
class TuneAD {
public:
  // Constructor
  TuneAD(const Rcpp::S4& LociMap,
         const Rcpp::S4& Pop,
         int nThreads_) : 
  nThreads(nThreads_) {
    // Assign class variables
    ploidy = Pop.slot("ploidy");
    nInd = Pop.slot("nInd");
    
    // Extract genotypes
    const arma::Col<int>& lociPerChr = LociMap.slot("lociPerChr");
    arma::uvec lociLoc = LociMap.slot("lociLoc");
    arma::Mat<unsigned char> genoMat = getGeno(
      Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
      lociPerChr, 
      lociLoc, 
      nThreads
    );
    
  }
  
private:
  arma::vec a, d, dd;
  double intercept, meanDD, varDD;
  arma::uword ploidy, nInd;
  int nThreads;
  
};

RCPP_MODULE(TuneAD_Module) {
  class_<TuneAD>("TuneAD")
  .constructor<const Rcpp::S4&, const Rcpp::S4&, int>()
  ;
}
