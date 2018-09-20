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
  if(hasGxe){
    output.set_size(2);
  }else{
    output.set_size(1);
  }
  arma::Mat<unsigned char> geno;
  geno = getGenoT(pop.slot("geno"), 
                  trait.slot("lociPerChr"),
                  trait.slot("lociLoc"));
  arma::vec a = trait.slot("addEff");
  double intercept = trait.slot("intercept");
  if(hasD){
    arma::vec d = trait.slot("domEff");
    output(0) = calcGvAD(geno, a, d, intercept, nThreads);
  }else{
    output(0) = calcGvA(geno, a, intercept, nThreads);
  }
  if(hasGxe){
    arma::vec g = trait.slot("gxeEff");
    double gxeInt = trait.slot("gxeInt");
    output(1) = calcGvA(geno, g, gxeInt, nThreads);
  }
  return output;
}

// A calculates breeding values and dominance deviations and genic
// variances. Additive and dominance genetic variances are calculated
// from breeding values and dominance deviations. Formulat accounts 
// for inbreeding in the population. Only works for ploidy=2.
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait, const Rcpp::S4& pop,
                        int nThreads){
  int nInd = pop.slot("nInd");
  arma::vec a = trait.slot("addEff");
  int nLoci = a.n_elem;
  arma::vec d(nLoci);
  if(trait.hasSlot("domEff")){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }else{
    d.zeros();
  }
  double intercept = trait.slot("intercept");
  arma::vec bv(nInd,arma::fill::zeros);
  arma::vec dd(nInd,arma::fill::zeros);
  arma::Mat<unsigned char> geno;
  geno = getGeno(pop.slot("geno"), 
                 trait.slot("lociPerChr"),
                 trait.slot("lociLoc"));
  arma::vec p(nLoci), q(nLoci), alpha(nLoci), mu(nLoci);
  arma::vec F(nLoci,arma::fill::zeros);
  arma::mat ddMat(3,nLoci);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(int i=0; i<nLoci; ++i){
    arma::vec genoFreq(3,arma::fill::zeros);
    for(int j=0; j<nInd; ++j){
      genoFreq(geno(j,i)) += 1;
    }
    p(i) = (genoFreq(2)+0.5*genoFreq(1))/accu(genoFreq);
    q(i) = 1-p(i);
    // 1-observed(het)/expect(het)
    if((p(i)>0.999999999) | (p(i)<0.000000001)){
      // Locus is fixed, no viable regression
      F(i) = 0;
      alpha(i) = 0;
      ddMat(0,i) = 0;
      ddMat(1,i) = 0;
      ddMat(2,i) = 0;
    }else{
      F(i) = 1-(genoFreq(1)/accu(genoFreq))/(2*p(i)*q(i));
      if(F(i)<-0.999999999){
        // Only heterozygotes, no viable regression
        alpha(i) = 0;
        ddMat(0,i) = 0;
        ddMat(1,i) = 0;
        ddMat(2,i) = 0;
      }else{
        double fFrac = (1-F(i))/(1+F(i));
        // a+d(q-p)(1-F)/(1+F)
        alpha(i) = a(i)+d(i)*(q(i)-p(i))*fFrac;
        // -2q(q+pF)(1-F)/(1+F)d
        ddMat(2,i) = -2*q(i)*(q(i)+p(i)*F(i))*fFrac*d(i);
        // (1-(p^2+q^2)(1-F)/(1+F)+2pqF(F-1)/(1+F))d
        ddMat(1,i) = (1-(p(i)*p(i)+q(i)*q(i))*fFrac+2*p(i)*q(i)*F(i)*(F(i)-1)/(1+F(i)))*d(i);
        // -2p(p+qF)(1-F)/(1+F)d
        ddMat(0,i) = -2*p(i)*(p(i)+q(i)*F(i))*fFrac*d(i);
      }
    }
    // 2pa+2pqd(1-F)
    mu(i) = 2*a(i)*p(i)+2*p(i)*q(i)*d(i)*(1-F(i));
  }
  arma::inplace_trans(geno);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(int i=0; i<nInd; ++i){
    for(int j=0; j<nLoci; ++j){
      bv(i) += (double(geno(j,i))-2*p(j))*alpha(j);
      dd(i) += ddMat(geno(j,i),j);
    }
  }
  intercept += accu(mu);
  // 2pq(1+F)alpha^2
  double genicVarA = 2.0*accu(p%q%arma::square(alpha)%(1+F));
  // 4pq(1-F)/(1+F)(p+Fq)(q+Fp)d^2
  double genicVarD = 4.0*accu(p%q%(1-F)/(1+F)%(p+F%q)%(q+F%p)%arma::square(d));
  return Rcpp::List::create(Rcpp::Named("bv")=bv,
                            Rcpp::Named("dd")=dd,
                            Rcpp::Named("genicVarA")=genicVarA,
                            Rcpp::Named("genicVarD")=genicVarD,
                            Rcpp::Named("mu")=intercept);
}
