// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates genetic values for traits with epistasis
arma::field<arma::vec> getGvE(const Rcpp::S4& trait, 
                              const Rcpp::S4& pop, 
                              int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::ivec& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
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
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc);
  
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
    for(arma::uword j=0; j<nInd; ++j){
      if(hasD){
        gv(j,tid) += a(E(i,0))*xa(genoMat(j,(E(i,0)))) + 
          d(E(i,0))*xd(genoMat(j,(E(i,0)))) + 
          a(E(i,1))*xa(genoMat(j,(E(i,1)))) + 
          d(E(i,1))*xd(genoMat(j,(E(i,1)))) + 
          E(i,2)*xa(genoMat(j,(E(i,0))))*xa(genoMat(j,(E(i,1))));
      }else{
        gv(j,tid) += a(E(i,0))*xa(genoMat(j,(E(i,0)))) + 
          a(E(i,1))*xa(genoMat(j,(E(i,1)))) + 
          E(i,2)*xa(genoMat(j,(E(i,0))))*xa(genoMat(j,(E(i,1))));
      }
      if(hasGxe){
        gxe(j,tid) += g(E(i,0))*xa(genoMat(j,(E(i,0)))) + 
          g(E(i,1))*xa(genoMat(j,(E(i,1))));
      }
    }
  }
  output(0) = sum(gv,1);
  if(hasGxe){
    output(1) = sum(gxe,1);
  }
  return output;
}

// Calculates genetic values for a trait
// Returns output in a list with length 1 or 2
//   The first item contains genetic values
//   The second item contains GxE effects (optional)
// [[Rcpp::export]]
arma::field<arma::vec> getGv(const Rcpp::S4& trait, 
                             const Rcpp::S4& pop, 
                             int nThreads){
  if(trait.hasSlot("epiEff")){
    return getGvE(trait, pop, nThreads);
  }
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nChr = pop.slot("nChr");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::field<arma::Cube<unsigned char> >& geno = pop.slot("geno");
  const arma::ivec& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  arma::mat gv(nInd,nChr),gxe;
  gv.fill(double(trait.slot("intercept"))/double(nChr));
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    output.set_size(2);
    output(0).set_size(nInd);
    output(1).set_size(nInd);
    gxe.set_size(nInd,nChr);
    gxe.fill(double(trait.slot("gxeInt"))/double(nChr));
  }else{
    output.set_size(1);
    output(0).set_size(nInd);
  }
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<nChr; ++i){
    if(lociPerChr(i)>0){ //Check for QTL
      arma::uword loc1=0,loc2;
      loc2 = arma::sum(lociPerChr(arma::span(0,i)))-1;
      if(i==0){
        loc1 = 0;
      }else{
        loc1 = arma::sum(lociPerChr(arma::span(0,i-1)));
      }
      arma::uvec chrLociLoc = lociLoc(arma::span(loc1,loc2))-1;
      arma::Mat<unsigned char> tmpGeno;
      tmpGeno = arma::sum(geno(i),1);
      tmpGeno = tmpGeno.rows(chrLociLoc).t();
      arma::vec eff(ploidy+1),gEff(ploidy+1);
      for(arma::uword j=loc1; j<(loc2+1); ++j){
        // Calculate genetic values
        eff = xa*a(j);
        if(hasD){
          eff += xd*d(j);
        }
        if(hasGxe){
          gEff = xa*g(j);
        }
        for(arma::uword k=0; k<nInd; ++k){
          gv(k,i) += eff(tmpGeno(k,j-loc1));
          if(hasGxe){
            gxe(k,i) += gEff(tmpGeno(k,j-loc1));
          }
        }
      }
    }
  }
  output(0) = sum(gv,1);
  if(hasGxe){
    output(1) = sum(gxe,1);
  }
  return output;
}
