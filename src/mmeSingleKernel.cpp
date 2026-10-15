#include "mme.h"

// Placeholders for the single trait, single kernel solvers, which are being
// reimplemented from published methods. They keep the exported signatures so
// that RcppExports does not change when the implementations arrive.

//' @title Solve RR-BLUP
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M){
  Rcpp::stop("solveRRBLUP: not yet reimplemented");
}

//' @title Solve Univariate Model
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Zu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K){
  Rcpp::stop("solveUVM: not yet reimplemented");
}
