#ifndef GETGV_H
#define GETGV_H

arma::vec calcGvA(const arma::Mat<unsigned char>& geno,
                  const arma::vec& a, double intercept);

arma::vec calcGvAD(const arma::Mat<unsigned char>& geno,
                   const arma::vec& a, const arma::vec& d,
                   double intercept);

#endif