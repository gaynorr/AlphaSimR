#ifndef GETGENO_H
#define GETGENO_H

arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                 const arma::Col<int>& lociPerChr,
                                 arma::uvec lociLoc, int nThreads);

arma::Mat<unsigned char> getMaternalGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads);

arma::Mat<unsigned char> getPaternalGeno(const arma::field<arma::Cube<unsigned char> >& geno, 
                                         const arma::Col<int>& lociPerChr,
                                         arma::uvec lociLoc, int nThreads);

arma::Mat<unsigned char> getOneHaplo(const arma::field<arma::Cube<unsigned char> >& geno, 
                                     const arma::Col<int>& lociPerChr,
                                     arma::uvec lociLoc, int haplo, int nThreads);

arma::mat genoToGenoA(const arma::Mat<unsigned char>& geno, 
                      arma::uword ploidy, int nThreads);

arma::mat genoToGenoD(const arma::Mat<unsigned char>& geno, 
                      arma::uword ploidy, int nThreads);

#endif