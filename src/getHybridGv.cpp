#include "alphasimr.h"

arma::field<arma::vec> getHybridGvE(const Rcpp::S4& trait, 
                                    const Rcpp::S4& females,
                                    arma::uvec& femaleParents,
                                    const Rcpp::S4& males,
                                    arma::uvec& maleParents,
                                    int nThreads){
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
  femaleParents -= 1; // R to C++
  maleParents -= 1; // R to C++
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = femaleParents.n_elem;
  arma::uword ploidyF = females.slot("ploidy");
  arma::uword ploidyM = males.slot("ploidy");
  // Calculate progeny ploidy without averaging
  // I'm not averaging because AlphaSimR scales for different ploidy
  arma::uword ploidy = ploidyF+ploidyM; 
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  arma::mat gv(nInd,nThreads),gxe;
  gv.fill(double(trait.slot("intercept"))/double(nThreads));
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    output.set_size(2);
    output(0).set_size(nInd);
    output(1).set_size(nInd);
    gxe.set_size(nInd,nThreads);
    gxe.fill(double(trait.slot("gxeInt"))/double(nThreads));
  }else{
    output.set_size(1);
    output(0).set_size(nInd);
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  arma::Mat<unsigned char> femaleGeno = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(females.slot("geno")), 
                                                lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> maleGeno = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(males.slot("geno")), 
                                              lociPerChr, lociLoc, nThreads);
  
  //Loop through loci pairs
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<E.n_rows; ++i){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    
    unsigned char geno1, geno2;
    for(arma::uword j=0; j<nInd; ++j){
      geno1 = femaleGeno(femaleParents(j),(E(i,0)))+maleGeno(maleParents(j),(E(i,0)));
      geno2 = femaleGeno(femaleParents(j),(E(i,1)))+maleGeno(maleParents(j),(E(i,1)));
      if(hasD){
        gv(j,tid) += a(E(i,0))*xa(geno1) + 
          d(E(i,0))*xd(geno1) + 
          a(E(i,1))*xa(geno2) + 
          d(E(i,1))*xd(geno2) + 
          E(i,2)*xa(geno1)*xa(geno2);
      }else{
        gv(j,tid) += a(E(i,0))*xa(geno1) + 
          a(E(i,1))*xa(geno2) + 
          E(i,2)*xa(geno1)*xa(geno2);
      }
      if(hasGxe){
        gxe(j,tid) += g(E(i,0))*xa(geno1) + 
          g(E(i,1))*xa(geno2);
      }
    }
  }
  output(0) = sum(gv,1);
  if(hasGxe){
    output(1) = sum(gxe,1);
  }
  return output;
  
}

// Calculates genetic values for cross between two inbreds
// Returns output in a list with length 1 or 2
//   The first item contains genetic values
//   The second item contains GxE effects (optional)
// [[Rcpp::export]]
arma::field<arma::vec> getHybridGv(const Rcpp::S4& trait, 
                                   const Rcpp::S4& females,
                                   arma::uvec femaleParents,
                                   const Rcpp::S4& males,
                                   arma::uvec maleParents,
                                   int nThreads){
  if(trait.hasSlot("epiEff")){
    return getHybridGvE(trait, females, femaleParents,
                        males, maleParents, nThreads);
  }
  femaleParents -= 1; // R to C++
  maleParents -= 1; // R to C++
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = femaleParents.n_elem;
  arma::uword ploidyF = females.slot("ploidy");
  arma::uword ploidyM = males.slot("ploidy");
  // Calculate progeny ploidy without averaging
  // I'm not averaging because AlphaSimR scales for different ploidy
  arma::uword ploidy = ploidyF+ploidyM; 
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  arma::mat gv(nInd,nThreads),gxe;
  gv.fill(double(trait.slot("intercept"))/double(nThreads));
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    output.set_size(2);
    output(0).set_size(nInd);
    output(1).set_size(nInd);
    gxe.set_size(nInd,nThreads);
    gxe.fill(double(trait.slot("gxeInt"))/double(nThreads));
  }else{
    output.set_size(1);
    output(0).set_size(nInd);
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  arma::Mat<unsigned char> femaleGeno = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(females.slot("geno")), 
                                                lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> maleGeno = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(males.slot("geno")), 
                                              lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<a.n_elem; ++i){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    arma::vec eff(ploidy+1),gEff(ploidy+1);
    eff = xa*a(i);
    if(hasD){
      eff += xd*d(i);
    }
    if(hasGxe){
      gEff = xa*g(i);
    }
    for(arma::uword j=0; j<nInd; ++j){
      gv(j,tid) += eff(femaleGeno(femaleParents(j),i)+
        maleGeno(maleParents(j),i));
      if(hasGxe){
        gxe(j,tid) += gEff(femaleGeno(femaleParents(j),i)+
          maleGeno(maleParents(j),i));
      }
    }
    
  }
  output(0) = sum(gv,1);
  if(hasGxe){
    output(1) = sum(gxe,1);
  }
  return output;
}