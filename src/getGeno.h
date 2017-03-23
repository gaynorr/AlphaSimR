#ifndef GETGENO_H
#define GETGENO_H

arma::Mat<unsigned char> getGeno(const Rcpp::S4& pop, const Rcpp::S4& lociMap);

arma::imat getDomGeno(const arma::Mat<unsigned char>& geno);

#endif