// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Calculates GEBVs for RRsol and RRDsol
// [[Rcpp::export]]
arma::mat gebvRR(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                 int nThreads){
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::mat a = sol.slot("markerEff");
  arma::mat fixEff = sol.slot("fixEff");
  arma::Mat<unsigned char> geno = getGeno(pop.slot("geno"), 
                                          sol.slot("lociPerChr"),
                                          sol.slot("lociLoc"),
                                          nThreads);
  arma::cube output(geno.n_rows,a.n_cols,nThreads);
  for(arma::uword i=0; i<a.n_cols; ++i){
    output.col(i).fill(fixEff(i)/double(nThreads));
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno.n_cols; ++j){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for(arma::uword i=0; i<geno.n_rows; ++i){
      output.slice(tid).row(i) += xa(geno(i,j))*a.row(j);
    }
  }
  return sum(output,2);
}

// Calculates GEGVs for RRDsol
// [[Rcpp::export]]
arma::mat gegvRRD(const Rcpp::S4& sol, const Rcpp::S4& pop,
                  int nThreads){
  arma::uword ploidy = pop.slot("ploidy");
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);
  arma::mat a = sol.slot("addEff");
  arma::mat d = sol.slot("domEff");
  arma::mat fixEff = sol.slot("fixEff");
  arma::Mat<unsigned char> geno = getGeno(pop.slot("geno"), 
                                          sol.slot("lociPerChr"),
                                          sol.slot("lociLoc"),
                                          nThreads);
  arma::cube output(geno.n_rows,a.n_cols,nThreads);
  for(arma::uword i=0; i<a.n_cols; ++i){
    output.col(i).fill(fixEff(i)/double(nThreads));
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno.n_cols; ++j){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for(arma::uword i=0; i<geno.n_rows; ++i){
      output.slice(tid).row(i) += xa(geno(i,j))*a.row(j) + xd(geno(i,j))*d.row(j);
    }
  }
  return sum(output,2);
}

// Assumes diploid
// [[Rcpp::export]]
arma::mat gebvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  bool female, int nThreads){
  arma::uword ploidy = 2;
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::mat fixEff = sol.slot("fixEff");
  arma::mat a;
  if(female){
    a = Rcpp::as<arma::mat>(sol.slot("femaleEff"));
  }else{
    a = Rcpp::as<arma::mat>(sol.slot("maleEff"));
  }
  arma::Mat<unsigned char> geno = getGeno(pop.slot("geno"), 
                                          sol.slot("lociPerChr"),
                                          sol.slot("lociLoc"),
                                          nThreads);
  arma::cube output(geno.n_rows,a.n_cols,nThreads);
  for(arma::uword i=0; i<a.n_cols; ++i){
    output.col(i).fill(fixEff(i)/double(nThreads));
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno.n_cols; ++j){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for(arma::uword i=0; i<geno.n_rows; ++i){
      output.slice(tid).row(i) += xa(geno(i,j))*a.row(j);
    }
  }
  return sum(output,2);
}

// Assumes diploid
// [[Rcpp::export]]
arma::mat gegvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  int nThreads){
  arma::uword ploidy = 1;
  double dP = double(ploidy);
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::mat fixEff = sol.slot("fixEff");
  arma::mat a1 = sol.slot("femaleEff");
  arma::mat a2 = sol.slot("maleEff");
  arma::Mat<unsigned char> geno1,geno2;
  geno1 = getOneHaplo(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       1,nThreads);
  geno2 = getOneHaplo(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       2,nThreads);
  arma::cube output(geno1.n_rows,a1.n_cols,nThreads);
  for(arma::uword i=0; i<a1.n_cols; ++i){
    output.col(i).fill(fixEff(i)/double(nThreads));
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno1.n_cols; ++j){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for(arma::uword i=0; i<geno1.n_rows; ++i){
      output.slice(tid).row(i) += xa(geno1(i,j))*a1.row(j) + 
        xa(geno2(i,j))*a2.row(j);
    }
  }
  return sum(output,2);
}

// Assumes diploid
// [[Rcpp::export]]
arma::mat gegvSCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  int nThreads){
  double dP = 1;
  arma::vec xa(2); // Genotype dossage
  for(arma::uword i=0; i<xa.n_elem; ++i)
    xa(i) = (double(i)-dP/2.0)*(2.0/dP);
  dP = 2;
  arma::vec xd(3);
  for(arma::uword i=0; i<xd.n_elem; ++i)
    xd(i) = double(i)*(dP-double(i))*(2.0/dP)*(2.0/dP);
  arma::mat fixEff = sol.slot("fixEff");
  arma::mat a1 = sol.slot("femaleEff");
  arma::mat a2 = sol.slot("maleEff");
  arma::mat d = sol.slot("d");
  arma::Mat<unsigned char> geno1,geno2,genoD;
  geno1 = getOneHaplo(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       1,nThreads);
  geno2 = getOneHaplo(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       2,nThreads);
  arma::cube output(geno1.n_rows,a1.n_cols,nThreads);
  for(arma::uword i=0; i<a1.n_cols; ++i){
    output.col(i).fill(fixEff(i)/double(nThreads));
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<geno1.n_cols; ++j){
    arma::uword tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    for(arma::uword i=0; i<geno1.n_rows; ++i){
      output.slice(tid).row(i) += xa(geno1(i,j))*a1.row(j) + 
        xa(geno2(i,j))*a2.row(j) + xd(geno1(i,j)+geno2(i,j))*d.row(j);
    }
  }
  return sum(output,2);
}
