// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates genetic values for a trait
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
  femaleParents -= 1; // R to C++
  maleParents -= 1; // R to C++
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nChr = females.slot("nChr");
  arma::uword nInd = femaleParents.n_elem;
  arma::uword ploidyF = females.slot("ploidy");
  arma::uword ploidyM = males.slot("ploidy");
  arma::uword ploidy = ploidyF+ploidyM;
  double dP = double(ploidy);
  const arma::field<arma::Cube<unsigned char> >& genoF = females.slot("geno");
  const arma::field<arma::Cube<unsigned char> >& genoM = males.slot("geno");
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
      arma::Mat<unsigned char> tmpGenoF, tmpGenoM;
      tmpGenoF = arma::sum(genoF(i),1);
      tmpGenoF = tmpGenoF.rows(chrLociLoc).t();
      tmpGenoM = arma::sum(genoM(i),1);
      tmpGenoM = tmpGenoM.rows(chrLociLoc).t();
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
          gv(k,i) += eff((tmpGenoF(femaleParents(k),j-loc1)+
            tmpGenoM(maleParents(k),j-loc1))/2);
          if(hasGxe){
            gxe(k,i) += gEff((tmpGenoF(femaleParents(k),j-loc1)+
              tmpGenoM(maleParents(k),j-loc1))/2);
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