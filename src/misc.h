#ifndef MISC_H
#define MISC_H

arma::uvec sampleInt(arma::uword n, arma::uword N);
arma::uword samplePoisson(double lambda);
double choose(double n, double k);
std::bitset<8> toBits(unsigned char byte);
unsigned char toByte(std::bitset<8> bits);

#endif
