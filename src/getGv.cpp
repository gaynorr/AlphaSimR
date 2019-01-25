// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates genetic values for a trait with only additive effects
arma::vec calcGvA(const arma::Mat<unsigned char>& geno,
                  const arma::vec& a, double intercept, 
                  int nThreads){
  arma::vec output(geno.n_cols);
  output.fill(intercept);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      output(i) += geno(j,i)*a(j);
    }
  }
  return output;
}

// Calculates genetic values for a trait with additive and dominance effects
arma::vec calcGvAD(const arma::Mat<unsigned char>& geno,
                   const arma::vec& a, const arma::vec& d,
                   double intercept, int nThreads){
  arma::vec output(geno.n_cols);
  output.fill(intercept);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      output(i) += geno(j,i)*a(j)+(1-abs(int(geno(j,i))-1))*d(j);
    }
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

double choose(double n, double k){ // n choose k
  if(k==0) return 1;
  return (n*choose(n-1,k-1))/k;
}

// A calculates breeding values and dominance deviations and genic
// variances. Additive and dominance genetic variances are calculated
// from breeding values and dominance deviations. Formulat accounts 
// for inbreeding in the population. Only works for ploidy=2.
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait, const Rcpp::S4& pop,
                        int nThreads){
  //Information from pop
  arma::uword nInd = pop.slot("nInd");
  arma::uword nChr  = pop.slot("nChr");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::field<arma::Cube<unsigned char> >& geno = pop.slot("geno");
  //Information from trait
  const arma::ivec& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::uword nLoci = a.n_elem;
  arma::vec d(nLoci);
  if(trait.hasSlot("domEff")){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }else{
    d.zeros();
  }
  double intercept = trait.slot("intercept");
  arma::mat bvMat(nInd,nChr,arma::fill::zeros); // "Breeding value"
  arma::mat ddMat(nInd,nChr,arma::fill::zeros); // Dominance deviation
  arma::mat gv_a(nInd,nChr,arma::fill::zeros); // Genetic value due to a
  arma::mat gv_d(nInd,nChr,arma::fill::zeros); // Genetic value due to d
  arma::vec genicA(nChr,arma::fill::zeros); // No LD
  arma::vec genicA2(nChr,arma::fill::zeros); // No LD and HWE
  arma::vec genicD(nChr,arma::fill::zeros); // No LD
  arma::vec genicD2(nChr,arma::fill::zeros); // No LD and HWE
  arma::vec mu(nChr,arma::fill::zeros); // Observed mean
  arma::vec eMu(nChr,arma::fill::zeros); // Expected mean with HWE
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
      arma::vec freq(ploidy+1), freqE(ploidy+1); // Genotype frequencies, observed and HWE
      arma::vec aEff(ploidy+1), dEff(ploidy+1), eff(ploidy+1); // Genetic values, additive and dominance
      arma::vec bv(ploidy+1), dd(ploidy+1), gv(ploidy+1); // Statistical values, additive and dominance
      arma::vec bvE(ploidy+1), ddE(ploidy+1); //Expected for random mating
      double gvMu, gvEMu, genoMu, p, q, dK, alpha, alphaE;
      for(arma::uword j=loc1; j<(loc2+1); ++j){
        // Calculate genotype and allele frequencies
        freq.zeros();
        for(arma::uword k=0; k<nInd; ++k){
          freq(tmpGeno(k,j-loc1)) += 1;
        }
        freq = freq/accu(freq);
        genoMu = accu(freq%x);
        p = genoMu/dP;
        q = 1-p;
        
        // Set effects, means and expected frequencies
        aEff = xa*a(j);
        dEff = xd*d(j);
        gv = aEff+dEff;
        gvMu = accu(freq%gv);
        mu(i) += gvMu;
        freqE.zeros();
        for(arma::uword k=0; k<(ploidy+1); ++k){
          dK = double(k);
          freqE(k) = choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
        }
        alpha = accu(freq%(gv-gvMu)%(x-genoMu))/
          accu(freq%(x-genoMu)%(x-genoMu));
        if(isinf(alpha)) alpha=0; //Check for divide by zero
        alphaE = accu(freqE%(gv-gvMu)%(x-genoMu))/
          accu(freqE%(x-genoMu)%(x-genoMu)); //Check for divide by zero
        if(isinf(alphaE)) alphaE=0;
        gvEMu =  accu(freqE%gv);
        eMu(i) += gvEMu;
        bv = (x-genoMu)*alpha; //Breeding values
        bvE = (x-genoMu)*alphaE; //Random mating breeding value
        dd = gv-bv-gvMu; //Dominance deviations (lack of fit)
        ddE = gv-bvE-gvEMu; //Random mating dominance deviation
        
        // Set genic variances
        genicA(i) += accu(freq%bv%bv);
        genicA2(i) += accu(freqE%bvE%bvE);
        genicD(i) += accu(freq%dd%dd);
        genicD2(i) += accu(freqE%ddE%ddE);
        
        // Set values for individuals
        for(arma::uword k=0; k<nInd; ++k){
          gv_a(k,i) += aEff(tmpGeno(k,j-loc1));
          gv_d(k,i) += dEff(tmpGeno(k,j-loc1));
          bvMat(k,i) += bv(tmpGeno(k,j-loc1));
          ddMat(k,i) += dd(tmpGeno(k,j-loc1));
        }
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("bv")=sum(bvMat,1),
                            Rcpp::Named("dd")=sum(ddMat,1),
                            Rcpp::Named("genicVarA")=accu(genicA),
                            Rcpp::Named("genicVarD")=accu(genicD),
                            Rcpp::Named("genicVarA2")=accu(genicA2),
                            Rcpp::Named("genicVarD2")=accu(genicD2),
                            Rcpp::Named("mu")=accu(mu)+intercept,
                            Rcpp::Named("inbreeding")=accu(eMu)-accu(mu),
                            Rcpp::Named("gv_a")=sum(gv_a,1),
                            Rcpp::Named("gv_d")=sum(gv_d,1),
                            Rcpp::Named("gv_mu")=intercept);
}
