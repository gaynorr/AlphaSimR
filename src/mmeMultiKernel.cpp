#include "mme.h"

// Placeholders for the single trait, multiple kernel solvers, which are being
// reimplemented from published methods. They keep the exported signatures so
// that RcppExports does not change when the implementations arrive.

//' @title Solve Multikernel RR-BLUP
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Mlist a list of M matrices
//' @param maxIter maximum number of iterations
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMK(arma::mat& y, arma::mat& X,
                         arma::field<arma::mat>& Mlist,
                         int maxIter=40){
  Rcpp::stop("solveRRBLUPMK: not yet reimplemented");
}

//' @title Solve Multikernel Model
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Zlist a list of Z matrices
//' @param Klist a list of K matrices
//' @param maxIter maximum number of iterations
//' @param tol tolerance for convergence
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMKM(arma::mat& y, arma::mat& X,
                    arma::field<arma::mat>& Zlist,
                    arma::field<arma::mat>& Klist,
                    int maxIter=40, double tol=1e-4){
  Rcpp::stop("solveMKM: not yet reimplemented");
}
