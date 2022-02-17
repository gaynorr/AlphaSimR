#ifndef MISC_H
#define MISC_H

arma::uword mapRow(const arma::uword& k, const arma::uword& n);
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n);
arma::uvec sampleInt(arma::uword n, arma::uword N);
arma::uword samplePoisson(double lambda);
arma::umat sampHalfDialComb(arma::uword nLevel, arma::uword n);
double choose(double n, double k);
std::bitset<8> toBits(unsigned char byte);
unsigned char toByte(std::bitset<8> bits);

#endif
