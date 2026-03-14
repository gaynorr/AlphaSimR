#ifndef ALPHASIMR_RNG_H
#define ALPHASIMR_RNG_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cstdint>
#include <dqrng_distribution.h>
#include <memory>

// AlphaSimR's RNG interface declarations live here.
// Non-template implementations live in rng.cpp.
namespace alphasimrRng {

using rngEngine = dqrng::random_64bit_generator;
using rngPtr = std::unique_ptr<rngEngine>;

dqrng::rng64_t createRng();
dqrng::rng64_t createRng(uint64_t seed);
rngPtr cloneStream(const dqrng::rng64_t &baseRng, uint64_t stream);

template <typename T> inline void shuffle(arma::Col<T> &x, rngEngine &rng) {
  std::shuffle(x.begin(), x.end(), rng);
}

double runif(rngEngine &rng);
arma::vec runifVec(arma::uword n, rngEngine &rng);
arma::vec rnormVec(arma::uword n, double mean, double sd, rngEngine &rng);
arma::vec rgammaVec(arma::uword n, double shape, double scale, rngEngine &rng);
arma::uvec sampleInt(arma::uword n, arma::uword N);
arma::uvec sampleInt(arma::uword n, arma::uword N, rngEngine &rng);
arma::uvec sampleInt(arma::uword n, arma::uword N, arma::uword seed);
arma::uword samplePoisson(double lambda, rngEngine &rng);

} // namespace alphasimrRng

#endif
