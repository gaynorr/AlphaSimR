#include "mme.h"

// Placeholders for the multiple trait, single kernel solvers, which are being
// reimplemented from published methods. They keep the exported signatures so
// that RcppExports does not change when the implementations arrive.

//' @title Solve Multivariate RR-BLUP
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Mu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param maxIter maximum number of iterations
//' @param tol tolerance for convergence
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X,
                         const arma::mat& M, int maxIter=1000,
                         double tol=1e-6){
  Rcpp::stop("solveRRBLUPMV: not yet reimplemented");
}

//' @title Solve Multivariate Model
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Zu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iterations
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K,
                    double tol=1e-6, int maxIter=1000){
  Rcpp::stop("solveMVM: not yet reimplemented");
}
