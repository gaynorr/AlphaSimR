#include "alphasimr.h"

/*
 * Class for tuning a trait with additive and dominance effects to achieve, 
 * as closely as possible, a desired level of dominance variance and 
 * inbreeding depression. The procedure should in the exact level of 
 * additive genetic variance and trait mean requested by the user. 
 */
class TuneAD {
public:
  
  // Constructor
  TuneAD(Rcpp::S4 LociMap,
         Rcpp::S4 Pop,
         double mean_,
         double varA_,
         double varD_,
         double inbrDepr_,
         int nThreads_) : 
  tarMean(mean_), tarVarA(varA_), 
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
      Rcpp::as<arma::field<arma::Cube<unsigned char> > >(Pop.slot("geno")), 
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
    genoFreq.zeros(ploidy+1, nLoci);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<nLoci; ++i){
      // Count genotypes
      for(arma::uword j=0; j<nInd; ++j){
        ++genoFreq(genoMat(j,i), i);
      }
      
      // Convert to frequency
      genoFreq.col(i) = genoFreq.col(i)/double(nInd);
      
      // Calculate genotype mean
      genoMu(i) = accu(genoFreq.col(i)%x);
      
      // Calculate inbreeding depression value
      double p = genoMu(i)/dP;
      double q = 1-p;
      // Expected heterozygosity at HWE
      // Not looping over first and last genotypes, because xd will be 0
      for(arma::uword k=1; k<(ploidy); ++k){
        double dK = double(k);
        hetHWE(i) += xd(k)*choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
      }
    }
  }
  
  // Objective function for optimization
  // Returns distance between desired inbrDepr and stdDevDom and 
  // calculated values using supplied meanDD and stdDevDD
  double objective(double meanDD_, double stdDevDD_){
    // Assign new values of meanDD and stdDevDD
    meanDD = meanDD_;
    stdDevDD = stdDevDD_;
    
    // Compute new values of d
    d = abs(a)%(domDegDev*stdDevDD + meanDD);
    
    // calcParam will scale additive and dominance effects
    // to achieve the desired varA
    calcParam();
    
    // Return distance between target and observed 
    // inbreeding depression and dominance variance (as standard deviation)
    return sqrt(std::pow(obsInbrDepr-tarInbrDepr, 2)+
                std::pow(sqrt(obsVarD)-sqrt(tarVarD), 2));
  }
  
  // Returns final values for a trait after setting meanDD and stdDevDD
  Rcpp::List finalize(double meanDD_, double stdDevDD_){
    // Assign new values of meanDD and stdDevDD
    meanDD = meanDD_;
    stdDevDD = stdDevDD_;
    
    // Compute new values of d
    d = abs(a)%(domDegDev*stdDevDD + meanDD);
    
    // calcParam will scale additive and dominance effects
    // to achieve the desired varA
    calcParam();
    
    // Return a's, d's and intercept of the trait
    // Report observed varD and inbreeding depression
    return Rcpp::List::create(Rcpp::Named("a")=a,
                              Rcpp::Named("d")=d,
                              Rcpp::Named("intercept")=intercept,
                              Rcpp::Named("varD")=obsVarD,
                              Rcpp::Named("inbrDepr")=obsInbrDepr);
  }
  
private:
  // Passed by constructor
  double tarMean, tarVarA, tarVarD, tarInbrDepr; // Targeted values
  int nThreads;
  
  // To be returned as output
  double intercept;
  arma::vec a;
  arma::vec d; // calculated from a and domDegDev
  double obsVarA, obsVarD, obsInbrDepr;
  
  // For internal calculations only
  double meanDD, stdDevDD;
  arma::vec x, xa, xd;
  arma::vec domDegDev;
  arma::vec genoMu, hetHWE;
  arma::Mat<unsigned char> genoMat;
  arma::mat genoFreq;
  arma::uword ploidy, nInd, nLoci;
  
  // Internal functions

  // Scales effects to hit target additive variance and records 
  // genetic parameters
  void calcParam(){
    // Allocate matrices for breeding values, dominance deviations, and means
    // Number of threads used for efficient parallel computing
    arma::mat bvMat(nInd, nThreads, arma::fill::zeros); // Breeding values
    arma::mat ddMat(nInd, nThreads, arma::fill::zeros); // Dominance deviations
    arma::vec muVec(nThreads, arma::fill::zeros); // Trait mean
    
    // Calculate breeding values and dominance deviations
    // Involves regressions for each locus
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<nLoci; ++i){
      
      // Assign thread ID
      arma::uword tid; 
#ifdef _OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      
      // Decompose genetic values into breeding values and dominance deviations
      arma::vec gv = xa*a(i) + xd*d(i);
      double gvMu = accu(genoFreq.col(i)%gv); // Mean genetic value
      muVec(tid) += gvMu; // Recording locus mean
      gv = gv - gvMu; // Centering genetic values
      arma::vec xc = x-genoMu(i); // Centered genotype dosage
      
      // Calculate average effect of an allele substitution (regression coefficient)
      double alpha = accu(genoFreq.col(i)%gv%xc)/
        accu(genoFreq.col(i)%xc%xc);
      
      // Calculate breeding values using alpha and dominance deviations using 
      // lack-of-fit
      arma::vec bv = xc*alpha;
      arma::vec dd = gv - bv;
      
      // Fill matrices for breeding values and dominance deviations
      // Accounts for the LD component of the variances
      for(arma::uword j=0; j<nInd; ++j){
        bvMat(j,tid) += bv(genoMat(j,i));
        ddMat(j,tid) += dd(genoMat(j,i));
      }
    }
    
    // Calculate additive and dominance genetic variances for population
    obsVarA = accu(sum(bvMat,1)%sum(bvMat,1));
    obsVarD = accu(sum(ddMat,1)%sum(ddMat,1));
    
    // Scale effects to hit target additive variance
    double scale = sqrt(tarVarA) / sqrt(obsVarA);
    a *= scale;
    d *= scale;
    obsVarA *= scale*scale;
    obsVarD *= scale*scale;
    
    // Calculate inbreeding depression and intercept for target mean
    obsInbrDepr = accu(hetHWE%d);
    intercept = tarMean - accu(muVec);
  }
};

RCPP_EXPOSED_CLASS(TuneAD)
RCPP_MODULE(TuneAD_module) {
  Rcpp::class_<TuneAD>("TuneAD")
  .constructor<Rcpp::S4, Rcpp::S4, double, double, double, double, int>()
  .method("objective", &TuneAD::objective)
  .method("finalize", &TuneAD::finalize)
  ;
}
