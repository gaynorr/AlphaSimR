// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"


// Retrieves GEBVs for RRsol
// [[Rcpp::export]]
arma::mat gebvRR(const Rcpp::S4& RRsol, const Rcpp::S4& pop, 
                 int nThreads){
  arma::mat a = RRsol.slot("markerEff");
  arma::Mat<unsigned char> geno;
  geno = getGenoT(pop.slot("geno"), 
                  RRsol.slot("lociPerChr"),
                  RRsol.slot("lociLoc"));
  arma::mat output(geno.n_cols,a.n_cols,arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      output.row(i) += geno(j,i)*a.row(j);
    }
  }
  return output;
}

// Retrieves GEGVs for RRDsol
// [[Rcpp::export]]
arma::mat gegvRRD(const Rcpp::S4& RRsol, const Rcpp::S4& pop,
                  int nThreads){
  arma::mat a = RRsol.slot("addEff");
  arma::mat d = RRsol.slot("domEff");
  double b = RRsol.slot("hetCov");
  arma::Mat<unsigned char> geno;
  geno = getGenoT(pop.slot("geno"), 
                  RRsol.slot("lociPerChr"),
                  RRsol.slot("lociLoc"));
  arma::mat output(geno.n_cols,a.n_cols,arma::fill::zeros);
  arma::vec het(geno.n_cols,arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      double dGeno = double(1-abs(int(geno(j,i))-1));
      output.row(i) += geno(j,i)*a.row(j)+dGeno*d.row(j);
      het(i) += dGeno;
    }
  }
  het = het/geno.n_rows;
  output += het*b;
  return output;
}


// [[Rcpp::export]]
arma::mat gebvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  bool female, int nThreads){
  arma::mat a;
  if(female){
    a = Rcpp::as<arma::mat>(sol.slot("femaleEff"));
  }else{
    a = Rcpp::as<arma::mat>(sol.slot("maleEff"));
  }
  arma::Mat<unsigned char> geno;
  geno = getGenoT(pop.slot("geno"), 
                  sol.slot("lociPerChr"),
                  sol.slot("lociLoc"));
  arma::mat output(geno.n_cols,a.n_cols,arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno.n_cols; ++i){
    for(arma::uword j=0; j<geno.n_rows; ++j){
      output.row(i) += geno(j,i)*a.row(j);
    }
  }
  return output;
}

// [[Rcpp::export]]
arma::mat gegvGCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  int nThreads){
  arma::mat a1 = sol.slot("femaleEff");
  arma::mat a2 = sol.slot("maleEff");
  arma::Mat<unsigned char> geno1,geno2;
  geno1 = getOneHaploT(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       1);
  geno2 = getOneHaploT(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       2);
  arma::mat output(geno1.n_cols,a1.n_cols,arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno1.n_cols; ++i){
    for(arma::uword j=0; j<geno1.n_rows; ++j){
      output.row(i) += geno1(j,i)*a1.row(j)+geno2(j,i)*a2.row(j);
    }
  }
  return 2*output;
}

// [[Rcpp::export]]
arma::mat gegvSCA(const Rcpp::S4& sol, const Rcpp::S4& pop, 
                  int nThreads){
  arma::mat a1 = sol.slot("a1");
  arma::mat a2 = sol.slot("a2");
  arma::mat d = sol.slot("d");
  double b = sol.slot("hetCov");
  arma::Mat<unsigned char> geno1,geno2,genoD;
  geno1 = getOneHaploT(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       1);
  geno2 = getOneHaploT(pop.slot("geno"), 
                       sol.slot("lociPerChr"),
                       sol.slot("lociLoc"), 
                       2);
  arma::mat output(geno1.n_cols,a1.n_cols,arma::fill::zeros);
  arma::vec het(geno1.n_cols,arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<geno1.n_cols; ++i){
    for(arma::uword j=0; j<geno1.n_rows; ++j){
      double dGeno = double(1-abs(int(geno1(j,i)+geno2(j,i))-1));
      output.row(i) += geno1(j,i)*a1.row(j)+geno2(j,i)*a2.row(j)+dGeno*d.row(j);
      het(i) += dGeno;
    }
  }
  het = het/geno1.n_rows;
  output += het*b;
  return output;
}
