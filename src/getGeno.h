#ifndef GETGENO_H
#define GETGENO_H

arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::ivec& lociPerChr,
                                 arma::uvec lociLoc);

arma::imat getDomGeno(const arma::Mat<unsigned char>& geno);

#endif