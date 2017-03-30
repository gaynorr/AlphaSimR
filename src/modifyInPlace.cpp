/*
 * The functions below are for modifying R objects in place.
 * 
 * The data type of oldValue must match exactly the data type in R or else
 * modifying in place will not occur.
 * 
 * For example, if assignMat is called on an integer matrix the function will
 * not modify in place. This is due to the function casting oldValue from 
 * integer matrix to a numeric matrix prior to assigning newValue.
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// [[Rcpp::export]]
void assignMat(arma::mat& oldValue, arma::mat& newValue){
  oldValue = newValue;
  return;
}

// [[Rcpp::export]]
void assignInt(arma::ivec& oldValue, arma::ivec& newValue){
  oldValue = newValue;
  return;
}

