#ifndef GETGENO_H
#define GETGENO_H

arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::ivec& lociPerChr,
                                 arma::uvec lociLoc);
arma::Mat<unsigned char> getGenoT(const arma::field<arma::Cube<unsigned char> >& geno, 
                                  const arma::ivec& lociPerChr,
                                  arma::uvec lociLoc);

arma::Mat<unsigned char> getOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                     const arma::ivec& lociPerChr,
                                     arma::uvec lociLoc, int haplo);

arma::Mat<unsigned char> getOneHaploT(const arma::field<arma::Cube<unsigned char> >& geno, 
                                      const arma::ivec& lociPerChr,
                                      arma::uvec lociLoc, int haplo);

#endif