#include "alphasimr.h"

// New class for tuning an AD trait
class TuneAD {
public:
  
  // Constructor
  TuneAD(const Rcpp::S4& LociMap,
         const Rcpp::S4& Pop,
         double mean_,
         double varA_,
         double varD_,
         double inbrDepr_,
         int nThreads_) : 
  tarMean(mean_), tarVarA(var_), 
  tarVarD(varD_), tarInbrDepr(inbrDepr_),
  nThreads(nThreads_) {
    
    // Assign class variables
    ploidy = Pop.slot("ploidy");
    double dP = double(ploidy);
    nInd = Pop.slot("nInd");
    x.set_size(ploidy+1);
    for(arma::uword i=0; i<x.n_elem; ++i)
      x(i) = double(i);
    xa = (x-dP/2.0)*(2.0/dP);
    xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
    
    // Extract genotypes
    const arma::Col<int>& lociPerChr = LociMap.slot("lociPerChr");
    nLoci = accu(lociPerChr);
    arma::uvec lociLoc = LociMap.slot("lociLoc");
    genoMat = getGeno(
      Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
      lociPerChr, 
      lociLoc, 
      nThreads
    );
    
    // Sample random deviates
    a.randn(nLoci);
    domDegDev.randn(nLoci);
    
    // Allocate vectors of length nLoci
    genoMu.set_size(nLoci);
    hetHWE.zeros(nLoci);
    
    // Calculate genotype frequencies
    genoFreq.zeros(ploidy+1, nLoci):
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<nLoci; ++i){
      // Count genotypes
      for(arma::uword j=0; j<nInd; ++j){
        ++genoFreq(genoMat(j,i), i);
      }
      
      // Convert to frequency
      genoFreq.col(i) = genFreq.col(i)/double(nInd);
      
      // Calculate genotype mean
      genoMu(i) = accu(genoFreq.col(i)%x);
      
      // Calculate inbreeding depression value
      double p = genoMu/dP;
      double q = 1-p;
      // Expected heterozygosity at HWE
      // Not looping over first and last genotypes, because xd will be 0
      for(arma::uword k=1; k<(ploidy); ++j){
        double dK = double(k);
        hetHWE(i) += xd(k)*choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
      }
    }
  }
  
  // Objective function for optimization
  // Returns distance between desired inbrDepr and stdDevDom and 
  // calculated values using supplied meanDD and stdDevDD
  double objective(meanDD_, stdDevDD_){
    meanDD = meanDD_;
    stdDevDD = stdDevDD_;
    d = fabs(a)*(domDegDev*stdDevDD + meanDD);
    calcParam();
    return sqrt(std::pow(obsInbrDepr-tarInbrDepr, 2)+
                std::pow(sqrt(obsVarD)-sqrt(tarVarD), 2));
  }
  
private:
  // Passed by constructor
  double tarMean, tarVarA, tarVarD, tarInbrDepr; // Targeted values
  int nThreads;
  
  // To be returned
  arma::vec a;
  arma::vec d; // calculated from a and domDegDev
  
  // For internal calculations
  double meanDD, stdDevDD;
  double obsVarA, obsVarD, obsInbrDepr;
  arma::vec x, xa, xd;
  arma::vec domDegDev;
  arma::vec genoMu, hetHWE;
  arma::Mat<unsigned char> genoMat;
  arma::mat genoFreq;
  arma::uword ploidy, nInd, nLoci;
  
  // Internal functions
  // Multistage function
  // Stage 1: calculate additive variance
  // Stage 2: scale additive effects for target variance
  // Stage 3: calculate variances and inbreeding depression
  void calcParam(){
    // Allocate matrices for breeding values and dominance deviations
    arma::mat bvMat(nInd, nThreads, arma::fill::zeros);
    arma::mat ddMat(nInd, nThreads, arma::fill::zeros);
    
    // Stage 1: Calculate breeding values with current a's and d's
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<a.n_elem; ++i){
      
      arma::uword tid; //Thread ID
#ifdef _OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      
      // Calculate breeding values
      arma::vec gv = xa*a(i) + xd*d(i);
      gv = gv - accu(genoFreq.col(i)%gv);
      double alpha = alpha = accu(genoFreq.col(i)%gv%(x-genoMu(i)))/
        accu(genoFreq.col(i)%(x-genoMu(i))%(x-genoMu(i)));
      arma::vec bv = (x-genoMu(i))*alpha;
      
      // Fill bvMat
      for(arma::uword j=0; j<nInd; ++j){
        bvMat(j,tid) += bv(genoMat(j,i));
      }
    }
    obsVarA = accu(sum(bvMat,1)%sum(bvMat,1));
    
    // Stage 2: rescale a's and d's
    a = a * sqrt(tarVarA) / sqrt(obsVarA);
    d = fabs(a) * (domDegDev*stdDevDD + meanDD);
    
    // Stage 3: calculate variances and inbreeding depression
    bvMat.zeros();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<a.n_elem; ++i){
      
      arma::uword tid; //Thread ID
#ifdef _OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      
      // Calculate breeding values
      arma::vec gv = xa*a(i) + xd*d(i);
      gv = gv - accu(genoFreq.col(i)%gv);
      double alpha = alpha = accu(genoFreq.col(i)%gv%(x-genoMu(i)))/
        accu(genoFreq.col(i)%(x-genoMu(i))%(x-genoMu(i)));
      arma::vec bv = (x-genoMu(i))*alpha;
      arma::vec dd = gv-bv;
      
      // Fill bvMat and ddMat
      for(arma::uword j=0; j<nInd; ++j){
        bvMat(j,tid) += bv(genoMat(j,i));
        ddMat(j,tid) += dd(genoMat(j,i));
      }
    }
    obsVarA = accu(sum(bvMat,1)%sum(bvMat,1));
    obsVarD = accu(sum(ddMat,1)%sum(ddMat,1));
    obsInbrDepr = accu(hetHWE%d);
  }
};

RCPP_MODULE(TuneAD_Module) {
  class_<TuneAD>("TuneAD")
  .constructor<const Rcpp::S4&, const Rcpp::S4&, double, double, double, double, int>()
  .method("objective", &TuneAD_Module::objective)
  ;
}
