#include <RcppArmadillo.h>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <cmath>

#include "rng.h"

// [[Rcpp::depends(RcppArmadillo)]]

// AlphaSimR's non-template RNG implementation lives here. See also rng.h.
//
// RNG policy in AlphaSimR:
// - R code uses R's RNG and can be reproduced via `set.seed()`.
// - MaCS is separate. It uses its own RNG internally, so R creates one integer
//   seed per chromosome, passes those seeds through the MaCS interface, and
//   AlphaSimR reuses the same seeds for post-MaCS site subsampling.
// - AlphaSimR's serial C++ RNG uses `createRng()` seeded from R's current RNG
//   state, so these serial C++ paths are still controlled by `set.seed()`.
// - OpenMP C++ code must not rely on a shared RNG. AlphaSimR's threaded C++
//   uses `createRng()` and `cloneStream()` to work across chromosome/work
//   items with stable, reproducible streams.
//
// AlphaSimR uses `dqrng` because it is fast and clonable, which makes it a
// good fit for both serial and OpenMP C++ code. C++ code could also use Rcpp,
// Armadillo, or R math helpers tied to R's RNG, but AlphaSimR avoids that for
// consistency and uses the RNG utilities implemented here.

namespace {

// Core implementation for sampling integers without replacement.
// n is the number of integers to return.
// N is the size of the source set, so sampled values range from 0 to N - 1.
// rng is a supplied RNG.
// Uses Jeffrey Scott Vitter's Method D.
arma::uvec sampleIntImpl(arma::uword n, arma::uword N,
                         alphasimrRng::rngEngine &rng) {
  arma::uvec output;
  if (n > N) {
    Rcpp::stop("sampleInt(): n must be <= N");
  }
  output.set_size(n);
  if (n == 0) {
    return output;
  }
  double q, v, x, y1, y2;
  arma::uword threshold = 13 * n;
  arma::uword S, limit, top, bottom;
  double u = alphasimrRng::runif(rng);
  v = std::exp(std::log(u) / double(n));
  q = double(N - n + 1);
  while ((n > 1) & (threshold < N)) {
    while (true) {
      while (true) {
        x = double(N) * (1 - v);
        S = std::floor(x);
        if (double(S) < q) {
          break;
        }
        u = alphasimrRng::runif(rng);
        v = std::exp(std::log(u) / double(n));
      }
      u = alphasimrRng::runif(rng);
      y1 = std::exp(std::log(u * double(N) / q) / double(n - 1));
      v = y1 * (1 - x / double(N)) * (q / (q - double(S)));
      if (v <= 1) {
        break;
      }
      y2 = 1;
      top = N - 1;
      if ((n - 1) > S) {
        bottom = N - n;
        limit = N - S;
      } else {
        bottom = N - S - 1;
        limit = N - n + 1;
      }
      for (arma::uword i = N - 1; i >= limit; --i) {
        y2 *= double(top) / double(bottom);
      }
      u = alphasimrRng::runif(rng);
      if ((double(N) / (double(N) - x)) >=
          (y1 * std::exp(std::log(y2) / double(n - 1)))) {
        v = std::exp(std::log(u) / double(n - 1));
        break;
      }
      v = std::exp(std::log(u) / double(n));
    }
    output(n - 1) = S + 1;
    N = N - S - 1;
    --n;
    q = double(N - n + 1);
    threshold -= 13;
  }
  if (n > 1) {
    top = N - n;
    while (n >= 2) {
      u = alphasimrRng::runif(rng);
      S = 0;
      q = double(top) / double(N);
      while (q > u) {
        ++S;
        --top;
        --N;
        q = (q * double(top)) / double(N);
      }
      output(n - 1) = S + 1;
      --N;
      --n;
    }
    u = alphasimrRng::runif(rng);
    output(0) = std::floor(u * N);
  } else {
    output(0) = std::floor(v * N);
  }
  return cumsum(output);
}

} // namespace

namespace alphasimrRng {

// Creates a clonable dqrng generator seeded from R's current RNG state.
// dqrng::generator<...>() draws the seed via get_seed_from_r(), so repeated
// runs with the same set.seed() and call sequence initialize the same dqrng
// state.
dqrng::rng64_t createRng() {
  return dqrng::generator<dqrng::xoshiro256plusplus>();
}

// Creates a clonable dqrng generator from an explicit integer seed.
// This is used when AlphaSimR needs RNG detached from R's current RNG state,
// such as seeded post-MaCS sampling inside an OpenMP region.
dqrng::rng64_t createRng(uint64_t seed) {
  return dqrng::generator<dqrng::xoshiro256plusplus>(seed);
}

// Clones a substream from an RNG.
// For xoshiro256++, clone(stream) uses a long-jump-derived substream,
// so callers should use stable logical work IDs such as chromosome indices
// rather than OpenMP thread IDs.
//
// Note that dqrng applies the jump as a loop, so this costs stream long
// jumps. Asking for a high numbered stream is therefore expensive, and
// building n streams by their ids costs O(n^2). Code that needs many
// streams should walk them instead, taking each one as clone(1) of the
// last, which costs one jump per stream and reaches the same states.
rngPtr cloneStream(const dqrng::rng64_t &baseRng, uint64_t stream) {
  return baseRng->clone(stream);
}

// Samples one U(0, 1) deviate from the supplied RNG.
double runif(rngEngine &rng) { return rng.uniform01(); }

// Samples n U(0, 1) deviates from the supplied RNG.
arma::vec runifVec(arma::uword n, rngEngine &rng) {
  arma::vec output(n);
  for (arma::uword i = 0; i < n; ++i) {
    output(i) = rng.uniform01();
  }
  return output;
}

// Samples n N(mean, sd^2) deviates from the supplied RNG.
// Use Boost here so normal and gamma draws share the same distribution layer.
arma::vec rnormVec(arma::uword n, double mean, double sd, rngEngine &rng) {
  arma::vec output(n);
  boost::random::normal_distribution<double> dist(mean, sd);
  for (arma::uword i = 0; i < n; ++i) {
    output(i) = dist(rng);
  }
  return output;
}

// Samples n Gamma(shape, scale) deviates from the supplied RNG.
// Use Boost here so normal and gamma draws share the same distribution layer.
arma::vec rgammaVec(arma::uword n, double shape, double scale, rngEngine &rng) {
  arma::vec output(n);
  boost::random::gamma_distribution<double> dist(shape, scale);
  for (arma::uword i = 0; i < n; ++i) {
    output(i) = dist(rng);
  }
  return output;
}

// Samples n integers from 0..N-1 without replacement, seeding from R's RNG.
// Safe for serial code; OpenMP code should use the rng or seed overloads.
arma::uvec sampleInt(arma::uword n, arma::uword N) {
  dqrng::rng64_t rng = createRng();
  return sampleIntImpl(n, N, *rng);
}

// Samples n integers from 0..N-1 without replacement using the supplied RNG.
arma::uvec sampleInt(arma::uword n, arma::uword N, rngEngine &rng) {
  return sampleIntImpl(n, N, rng);
}

// Samples n integers from 0..N-1 without replacement using the supplied seed.
// This is used primarily in the MaCS path.
arma::uvec sampleInt(arma::uword n, arma::uword N, arma::uword seed) {
  dqrng::rng64_t rng = createRng(seed);
  return sampleIntImpl(n, N, *rng);
}

// Samples one Poisson(lambda) deviate using the supplied RNG.
// Uses Knuth's method, which is suitable for the small lambda values used in
// AlphaSimR.
arma::uword samplePoisson(double lambda, rngEngine &rng) {
  double p = 1;
  double L = std::exp(-lambda);
  arma::uword k = 0;
  do {
    ++k;
    p *= runif(rng);
  } while (p > L);
  return k - 1;
}

} // namespace alphasimrRng

// [[Rcpp::export]]
arma::imat rngDiagnosticsSampleInt(arma::uword n, arma::uword N,
                                   arma::uword reps, uint64_t seed) {
  dqrng::rng64_t rng = alphasimrRng::createRng(seed);
  arma::imat output(reps, n);

  for (arma::uword rep = 0; rep < reps; ++rep) {
    arma::uvec draw = alphasimrRng::sampleInt(n, N, *rng);
    output.row(rep) = arma::conv_to<arma::irowvec>::from(draw.t());
  }

  return output;
}

// [[Rcpp::export]]
arma::ivec rngDiagnosticsSamplePoisson(double lambda, arma::uword reps,
                                       uint64_t seed) {
  dqrng::rng64_t rng = alphasimrRng::createRng(seed);
  arma::ivec output(reps);

  for (arma::uword rep = 0; rep < reps; ++rep) {
    output(rep) = static_cast<int>(alphasimrRng::samplePoisson(lambda, *rng));
  }

  return output;
}
