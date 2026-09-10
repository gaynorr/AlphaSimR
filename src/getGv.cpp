#include "alphasimr.h"

// Calculates genetic values for genomic predictions using parental origin
arma::field<arma::vec> getGvA2(const Rcpp::S4& trait, 
                               const Rcpp::S4& pop, 
                               int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a1,a2,d;
  a1 = Rcpp::as<arma::vec>(trait.slot("addEff"));
  a2 = Rcpp::as<arma::vec>(trait.slot("addEffMale"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(a1.n_elem);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks);
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  output.set_size(1);
  output(0).set_size(nInd);
  // Half ploidy for xa
  arma::vec xa(ploidy/2+1);
  for(arma::uword i=0; i<xa.n_elem; ++i)
    xa(i) = (double(i)-dP/4.0)*(4.0/dP);
  // Full ploidy level for xd
  arma::vec xd(ploidy+1);
  for(arma::uword i=0; i<xd.n_elem; ++i)
    xd(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
  
  
  arma::Mat<unsigned char> maternalGeno = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> paternalGeno = getPaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(a1.n_elem, nBlocks, tid);
    arma::uword itemEnd = blockStart(a1.n_elem, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      arma::vec aEff1,aEff2,dEff;
      aEff1 = xa*a1(i);
      aEff2 = xa*a2(i);
      if(hasD){
        dEff = xd*d(i);
      }
      for(arma::uword j=0; j<nInd; ++j){
        gv(j,tid) += aEff1(maternalGeno(j,i)) + aEff2(paternalGeno(j,i));
        if(hasD){
          gv(j,tid) += dEff(maternalGeno(j,i)+paternalGeno(j,i));
        }
      }
    }
  }
  output(0) = sum(gv,1);
  return output;
}

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
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
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
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(E.n_rows);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks),gxe;
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    output.set_size(2);
    output(0).set_size(nInd);
    output(1).set_size(nInd);
    gxe.set_size(nInd,nBlocks);
    gxe.fill(double(trait.slot("gxeInt"))/double(nBlocks));
  }else{
    output.set_size(1);
    output(0).set_size(nInd);
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  
  //Loop through loci pairs
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(E.n_rows, nBlocks, tid);
    arma::uword itemEnd = blockStart(E.n_rows, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
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
  if(trait.hasSlot("addEffMale")){
    // Genomic prediction
    return getGvA2(trait, pop, nThreads);
  }
  if(trait.hasSlot("epiEff")){
    return getGvE(trait, pop, nThreads);
  }
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(a.n_elem);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks),gxe;
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    output.set_size(2);
    output(0).set_size(nInd);
    output(1).set_size(nInd);
    gxe.set_size(nInd,nBlocks);
    gxe.fill(double(trait.slot("gxeInt"))/double(nBlocks));
  }else{
    output.set_size(1);
    output(0).set_size(nInd);
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  arma::Mat<unsigned char> genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                             lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(a.n_elem, nBlocks, tid);
    arma::uword itemEnd = blockStart(a.n_elem, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      arma::vec eff(ploidy+1),gEff(ploidy+1);
      eff = xa*a(i);
      if(hasD){
        eff += xd*d(i);
      }
      if(hasGxe){
        gEff = xa*g(i);
      }
      for(arma::uword j=0; j<nInd; ++j){
        gv(j,tid) += eff(genoMat(j,i));
        if(hasGxe){
          gxe(j,tid) += gEff(genoMat(j,i));
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

// Basic implementation of getGv for multiple traits
arma::field<arma::vec> getGvIndexStd(const Rcpp::S4& trait, 
                                     arma::Mat<unsigned char>& genoMat, 
                                     arma::Col<int> qtlIndex,
                                     arma::uword ploidy,
                                     int nThreads){
  qtlIndex = qtlIndex-1; // R to C++
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = genoMat.n_rows;
  arma::uword nGeno = ploidy+1;
  double dP = double(ploidy);
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  arma::uword nLoci = a.n_elem;
  arma::vec x(nGeno); // Genotype dosage
  for(arma::uword i=0; i<nGeno; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  // What a locus contributes depends only on its dosage, so the value of
  // each dosage is worked out once per locus and looked up afterwards
  arma::mat effTab(nGeno,nLoci);
  for(arma::uword i=0; i<nLoci; ++i){
    if(hasD){
      effTab.col(i) = xa*a(i)+xd*d(i);
    }else{
      effTab.col(i) = xa*a(i);
    }
  }
  arma::mat gxeTab;
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
    gxeTab.set_size(nGeno,nLoci);
    for(arma::uword i=0; i<nLoci; ++i){
      gxeTab.col(i) = xa*g(i);
    }
  }
  
  // Written into directly, so no accumulator is needed and nothing has to be
  // summed at the end
  if(hasGxe){
    output.set_size(2);
  }else{
    output.set_size(1);
  }
  output(0).set_size(nInd);
  output(0).fill(double(trait.slot("intercept")));
  double* gvPtr = output(0).memptr();
  double* gxePtr = NULL;
  if(hasGxe){
    output(1).set_size(nInd);
    output(1).fill(double(trait.slot("gxeInt")));
    gxePtr = output(1).memptr();
  }
  
  // Work is split into blocks of individuals. Each block owns its own part
  // of the output, so the split does not depend on the thread count and an
  // individual always sums its loci in the same order. A block's genetic
  // values also stay in cache while the loci are swept, instead of the whole
  // output being walked once per locus.
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword i=0; i<nLoci; ++i){
      const unsigned char* geno = genoMat.colptr(arma::uword(qtlIndex(i)));
      const double* eff = effTab.colptr(i);
      for(arma::uword j=indStart; j<indEnd; ++j){
        gvPtr[j] += eff[geno[j]];
      }
      if(hasGxe){
        const double* gEff = gxeTab.colptr(i);
        for(arma::uword j=indStart; j<indEnd; ++j){
          gxePtr[j] += gEff[geno[j]];
        }
      }
    }
  }
  return output;
}

// Calculates genetic values for multiple traits with epistasis
arma::field<arma::vec> getGvIndexE(const Rcpp::S4& trait, 
                                   arma::Mat<unsigned char>& genoMat, 
                                   arma::Col<int> qtlIndex,
                                   arma::uword ploidy,
                                   int nThreads){
  qtlIndex = qtlIndex-1; // R to C++
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  bool hasGxe = trait.hasSlot("gxeEff");
  arma::uword nInd = genoMat.n_rows;
  arma::uword nGeno = ploidy+1;
  double dP = double(ploidy);
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
  arma::vec a,d,g;
  a = Rcpp::as<arma::vec>(trait.slot("addEff"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  if(hasGxe){
    g = Rcpp::as<arma::vec>(trait.slot("gxeEff"));
  }
  arma::uword nPair = E.n_rows;
  arma::vec x(nGeno); // Genotype dosage
  for(arma::uword i=0; i<nGeno; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  
  // Every term of the sum depends only on the two dosages, so the value of
  // each pair of dosages is worked out once per locus pair and looked up
  // afterwards. The columns holding the two loci are recorded at the same
  // time, so neither the epistasis matrix nor the QTL index is read again.
  arma::mat effTab(nGeno*nGeno,nPair), gxeTab;
  if(hasGxe){
    gxeTab.set_size(nGeno*nGeno,nPair);
  }
  arma::uvec col1(nPair), col2(nPair);
  for(arma::uword i=0; i<nPair; ++i){
    arma::uword l1 = arma::uword(E(i,0));
    arma::uword l2 = arma::uword(E(i,1));
    col1(i) = arma::uword(qtlIndex(l1));
    col2(i) = arma::uword(qtlIndex(l2));
    double a1 = a(l1), a2 = a(l2), e12 = E(i,2);
    double d1 = 0.0, d2 = 0.0, g1 = 0.0, g2 = 0.0;
    if(hasD){
      d1 = d(l1);
      d2 = d(l2);
    }
    if(hasGxe){
      g1 = g(l1);
      g2 = g(l2);
    }
    double* eff = effTab.colptr(i);
    double* gEff = NULL;
    if(hasGxe){
      gEff = gxeTab.colptr(i);
    }
    for(arma::uword u=0; u<nGeno; ++u){
      for(arma::uword v=0; v<nGeno; ++v){
        if(hasD){
          eff[u+v*nGeno] = a1*xa(u)+d1*xd(u)+a2*xa(v)+d2*xd(v)+e12*xa(u)*xa(v);
        }else{
          eff[u+v*nGeno] = a1*xa(u)+a2*xa(v)+e12*xa(u)*xa(v);
        }
        if(hasGxe){
          gEff[u+v*nGeno] = g1*xa(u)+g2*xa(v);
        }
      }
    }
  }
  
  // Written into directly, so no accumulator is needed and nothing has to be
  // summed at the end
  if(hasGxe){
    output.set_size(2);
  }else{
    output.set_size(1);
  }
  output(0).set_size(nInd);
  output(0).fill(double(trait.slot("intercept")));
  double* gvPtr = output(0).memptr();
  double* gxePtr = NULL;
  if(hasGxe){
    output(1).set_size(nInd);
    output(1).fill(double(trait.slot("gxeInt")));
    gxePtr = output(1).memptr();
  }
  
  // Work is split into blocks of individuals, as in getGvIndexStd
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  //Loop through loci pairs
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword i=0; i<nPair; ++i){
      const unsigned char* geno1 = genoMat.colptr(col1(i));
      const unsigned char* geno2 = genoMat.colptr(col2(i));
      const double* eff = effTab.colptr(i);
      for(arma::uword j=indStart; j<indEnd; ++j){
        gvPtr[j] += eff[geno1[j]+geno2[j]*nGeno];
      }
      if(hasGxe){
        const double* gEff = gxeTab.colptr(i);
        for(arma::uword j=indStart; j<indEnd; ++j){
          gxePtr[j] += gEff[geno1[j]+geno2[j]*nGeno];
        }
      }
    }
  }
  return output;
}

// Calculates genetic values for genomic predictions using parental origin
arma::field<arma::vec> getGvIndexA2(const Rcpp::S4& trait, 
                                    const Rcpp::S4& pop, 
                                    int nThreads){
  arma::field<arma::vec> output;
  bool hasD = trait.hasSlot("domEff");
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a1,a2,d;
  a1 = Rcpp::as<arma::vec>(trait.slot("addEff"));
  a2 = Rcpp::as<arma::vec>(trait.slot("addEffMale"));
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  // Accumulators are indexed by work block, not by thread, so that
  // their number and the order they are summed in do not depend on
  // how many threads are available
  arma::uword nBlocks = countBlocks(a1.n_elem);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::mat gv(nInd,nBlocks);
  gv.fill(double(trait.slot("intercept"))/double(nBlocks));
  output.set_size(1);
  output(0).set_size(nInd);
  // Half ploidy for xa
  arma::vec xa(ploidy/2+1);
  for(arma::uword i=0; i<xa.n_elem; ++i)
    xa(i) = (double(i)-dP/4.0)*(4.0/dP);
  // Full ploidy level for xd
  arma::vec xd(ploidy+1);
  for(arma::uword i=0; i<xd.n_elem; ++i)
    xd(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
  
  
  arma::Mat<unsigned char> maternalGeno = getMaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  arma::Mat<unsigned char> paternalGeno = getPaternalGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                                                          lociPerChr, lociLoc, nThreads);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(a1.n_elem, nBlocks, tid);
    arma::uword itemEnd = blockStart(a1.n_elem, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      arma::vec aEff1,aEff2,dEff;
      aEff1 = xa*a1(i);
      aEff2 = xa*a2(i);
      if(hasD){
        dEff = xd*d(i);
      }
      for(arma::uword j=0; j<nInd; ++j){
        gv(j,tid) += aEff1(maternalGeno(j,i)) + aEff2(paternalGeno(j,i));
        if(hasD){
          gv(j,tid) += dEff(maternalGeno(j,i)+paternalGeno(j,i));
        }
      }
    }
  }
  output(0) = sum(gv,1);
  return output;
}

// [[Rcpp::export]]
Rcpp::List getGvIndex(const Rcpp::S4& pop, 
                      const Rcpp::List& traitList, 
                      const Rcpp::S4& activeQtl,
                      const arma::field<arma::Col<int> >& qtlIndex,
                      arma::uword nTraits,
                      int nThreads){
  
  // Get summary variables
  arma::uword nInd = pop.slot("nInd");
  arma::uword ploidy = pop.slot("ploidy");
  
  // Get genotypes
  arma::Mat<unsigned char> genoMat; 
  genoMat = getGeno(Rcpp::as<arma::field<arma::Cube<unsigned char> > >(pop.slot("geno")), 
                    Rcpp::as<arma::Col<int> >(activeQtl.slot("lociPerChr")),
                    Rcpp::as<arma::uvec >(activeQtl.slot("lociLoc")),
                    nThreads);
  
  // Format output
  arma::mat gv(nInd, nTraits);
  arma::field<arma::vec> gxe(nTraits);
  
  
  // Loop over traits
  for(arma::uword i=0; i<nTraits; ++i){
    const Rcpp::S4& trait = traitList[i];
    // An arma field, not an Rcpp list, so the values do not make a round
    // trip through R objects on the way back
    arma::field<arma::vec> output;
    if(trait.hasSlot("addEffMale")){
      // Trait with separate male and female effects (from genomic prediction)
      // Unlikely to be use
      output = getGvIndexA2(trait, pop, nThreads);
    }else if(trait.hasSlot("epiEff")){
      // Trait with epistasis
      output = getGvIndexE(trait, genoMat, qtlIndex(i), ploidy, nThreads);
    }else{
      // All other traits
      output = getGvIndexStd(trait, genoMat, qtlIndex(i), ploidy, nThreads);
    }
    
    gv.col(i) = output(0);
    
    if(output.n_elem == 2){
      // Has GxE
      gxe(i) = output(1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("gv") = gv, 
                            Rcpp::Named("gxe") = gxe);
}