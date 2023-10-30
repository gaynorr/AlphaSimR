#include "alphasimr.h"

// Sets up the list of arguments needed for optimization
// [[Rcpp::export]]
Rcpp::List argAltAD(Rcpp::S4 LociMap,
                    Rcpp::S4 Pop,
                    double mean,
                    double varA,
                    double varD,
                    double inbrDepr,
                    int nThreads){
  
  // Create ploidy specific genotype dosage variables
  arma::uword ploidy = Pop.slot("ploidy");
  double dP = double(ploidy);
  arma::uword nInd = Pop.slot("nInd");
  arma::vec x(ploidy+1);
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  // Extract loci information and genotypes
  const arma::Col<int>& lociPerChr = LociMap.slot("lociPerChr");
  arma::uword nLoci = accu(lociPerChr);
  arma::uvec lociLoc = LociMap.slot("lociLoc");
  arma::Mat<unsigned char> genoMat = getGeno(
    Rcpp::as<arma::field<arma::Cube<unsigned char> > >(Pop.slot("geno")), 
    lociPerChr, 
    lociLoc, 
    nThreads
  );
  
  // Calculate genotype frequencies
  arma::mat genoFreq(ploidy+1, nLoci, arma::fill::zeros);
  arma::vec genoMu(nLoci), hetHWE(nLoci);
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
  
  // Sample random deviates
  arma::vec a(nLoci, arma::fill::randn);
  arma::vec domDegDev(nLoci, arma::fill::randn);
  
  return Rcpp::List::create(Rcpp::Named("x")=x,
                            Rcpp::Named("xa")=xa,
                            Rcpp::Named("xd")=xd,
                            Rcpp::Named("genoMat")=genoMat,
                            Rcpp::Named("genoFreq")=genoFreq,
                            Rcpp::Named("genoMu")=genoMu,
                            Rcpp::Named("hetHWE")=hetHWE,
                            Rcpp::Named("a")=a,
                            Rcpp::Named("domDegDev")=domDegDev,
                            Rcpp::Named("mean")=mean,
                            Rcpp::Named("varA")=varA,
                            Rcpp::Named("varD")=varD,
                            Rcpp::Named("inbrDepr")=inbrDepr,
                            Rcpp::Named("nThreads")=nThreads);
}

// The objective function for optimization
// [[Rcpp::export]]
double objAltAD(arma::vec input, const Rcpp::List& args){
  // Assign new values of meanDD and stdDevDD
  double meanDD = input(0);
  double stdDevDD = input(1);
  
  // Access variables from args
  const arma::Mat<unsigned char>& genoMat = args["genoMat"];
  arma::uword nInd = genoMat.n_rows;
  arma::uword nLoci = genoMat.n_cols;
  const arma::mat& genoFreq = args["genoFreq"];
  const arma::vec& a = args["a"];
  const arma::vec& domDegDev = args["domDegDev"];
  const arma::vec& x = args["x"];
  const arma::vec& xa = args["xa"];
  const arma::vec& xd = args["xd"];
  const arma::vec& genoMu = args["genoMu"];
  const arma::vec& hetHWE = args["hetHWE"];
  double varA = args["varA"];
  double varD = args["varD"];
  double inbrDepr = args["inbrDepr"];
  int nThreads = args["nThreads"];
  
  // Calculate d
  arma::vec d = abs(a)%(domDegDev*stdDevDD + meanDD);
  
  // Allocate matrices for breeding values, dominance deviations, and means
  // Number of threads used for efficient parallel computing
  arma::mat bvMat(nInd, nThreads, arma::fill::zeros); // Breeding values
  arma::mat ddMat(nInd, nThreads, arma::fill::zeros); // Dominance deviations
  
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
  double obsVarA = accu(sum(bvMat,1)%sum(bvMat,1)) / nInd;
  double obsVarD = accu(sum(ddMat,1)%sum(ddMat,1)) / nInd;
  
  // Scale effects to hit target additive variance
  double scale = sqrt(varA) / sqrt(obsVarA);
  d *= scale;
  obsVarD *= scale*scale;
  
  // Calculate inbreeding depression and intercept for target mean
  double obsInbrDepr = accu(hetHWE%d);
  
  // Return distance between target and observed 
  // inbreeding depression and dominance variance (as standard deviation)
  return sqrt(std::pow(obsInbrDepr-inbrDepr, 2) +
              std::pow(sqrt(obsVarD)-sqrt(varD), 2));
}

// Calculates the a and d effects and the intercept
// [[Rcpp::export]]
Rcpp::List finAltAD(arma::vec input, const Rcpp::List& args){
  // Assign new values of meanDD and stdDevDD
  double meanDD = input(0);
  double stdDevDD = input(1);
  
  // Access variables from args
  const arma::Mat<unsigned char>& genoMat = args["genoMat"];
  arma::uword nInd = genoMat.n_rows;
  arma::uword nLoci = genoMat.n_cols;
  const arma::mat& genoFreq = args["genoFreq"];
  arma::vec a = args["a"];
  const arma::vec& domDegDev = args["domDegDev"];
  const arma::vec& x = args["x"];
  const arma::vec& xa = args["xa"];
  const arma::vec& xd = args["xd"];
  const arma::vec& genoMu = args["genoMu"];
  const arma::vec& hetHWE = args["hetHWE"];
  double varA = args["varA"];
  double mean = args["mean"];
  int nThreads = args["nThreads"];
  
  // Calculate d
  arma::vec d = abs(a)%(domDegDev*stdDevDD + meanDD);
  
  // Allocate matrices for breeding values, dominance deviations, and means
  // Number of threads used for efficient parallel computing
  arma::mat bvMat(nInd, nThreads, arma::fill::zeros); // Breeding values
  arma::mat ddMat(nInd, nThreads, arma::fill::zeros); // Dominance deviations
  
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
  double obsVarA = accu(sum(bvMat,1)%sum(bvMat,1)) / nInd;
  double obsVarD = accu(sum(ddMat,1)%sum(ddMat,1)) / nInd;
  
  // Scale effects to hit target additive variance
  double scale = sqrt(varA) / sqrt(obsVarA);
  a *= scale;
  d *= scale;
  obsVarD *= scale*scale;
  
  // Calculate inbreeding depression and intercept for target mean
  double obsInbrDepr = accu(hetHWE%d);
  
  // Calculate GV
  arma::mat gvMat(nInd, nThreads, arma::fill::zeros); // Genetic values
  
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
    
    // Fill matrices for breeding values and dominance deviations
    // Accounts for the LD component of the variances
    for(arma::uword j=0; j<nInd; ++j){
      gvMat(j,tid) += gv(genoMat(j,i));
    }
  }
  
  gvMat = sum(gvMat,1);
  
  double intercept = mean - accu(gvMat) / nInd;
  
  return Rcpp::List::create(Rcpp::Named("a")=a,
                            Rcpp::Named("d")=d,
                            Rcpp::Named("intercept")=intercept,
                            Rcpp::Named("meanDD")=meanDD,
                            Rcpp::Named("varDD")=stdDevDD*stdDevDD,
                            Rcpp::Named("inbrDepr")=obsInbrDepr,
                            Rcpp::Named("varD")=obsVarD,
                            Rcpp::Named("varG")=accu(arma::cov(gvMat,1)));
}
