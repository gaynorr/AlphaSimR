#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
/**
 * First and raw second moments of total genealogical branch length.
 *
 * Entries are indexed by the starting number of ancestral lineages. Index 0
 * is unused, and index 1 stays zero because the genealogy ends at its most
 * recent common ancestor. Means use units of 2*N0 generations; second moments
 * use the square of those units. The second vector stores E[T^2], not Var(T).
 */
struct TreeMoments {
  std::vector<double> mean;
  std::vector<double> second;
  /** Allocate zero moments for lineage counts 0 through n. */
  explicit TreeMoments(int n) : mean(n + 1, 0.0), second(n + 1, 0.0) {}
};

/**
 * Advance the backward moment equations by one implicit Euler step.
 *
 * Extending an epoch toward the present adds branch length at rate k when k
 * lineages remain. With q_k = k*(k-1)/(2*size), the moment equations are
 *   dM_k/ds = k + q_k*(M_(k-1) - M_k),
 *   dQ_k/ds = 2*k*M_k + q_k*(Q_(k-1) - Q_k).
 * Here s measures how far the epoch has been extended toward the present,
 * M_k = E[T | k], and Q_k = E[T^2 | k]. Propagating the raw second moment
 * includes correlations between the times spent at different lineage counts.
 *
 * The implicit system is triangular: ascending lineage counts allow each
 * updated state to use the already updated state k-1. This handles the large
 * coalescence rates at high k without constructing a dense transition matrix.
 *
 * @param value Moments at the older end of the step, including all older epochs.
 * @param step Positive step length in units of 2*N0 generations.
 * @param size Positive population size relative to N0, constant over the step.
 * @return Updated first and raw second moments for all lineage counts.
 *         Runtime and temporary storage are both O(n).
 */
TreeMoments implicitStep(const TreeMoments &value, double step, double size) {
  int n = value.mean.size() - 1;
  TreeMoments next(n);
  for (int k = 2; k <= n; ++k) {
    double rate = k * (k - 1.0) / (2.0 * size);
    double weight = step * rate;
    next.mean[k] =
        (value.mean[k] + step * k + weight * next.mean[k - 1]) / (1.0 + weight);
    next.second[k] = (value.second[k] + 2.0 * step * k * next.mean[k] +
                      weight * next.second[k - 1]) /
                     (1.0 + weight);
  }
  return next;
}

/**
 * Obtain a second-order step by extrapolating implicit Euler solutions.
 *
 * Combining two half steps and one full step as 2*half - full cancels the
 * leading error term. This function does not choose or accept a step size;
 * the outer solver compares full and half second-order steps for that purpose.
 *
 * @param value Moments at the older end of the step.
 * @param step Positive step length in units of 2*N0 generations.
 * @param size Positive population size relative to N0 within the step.
 * @return Extrapolated first and raw second moments for all lineage counts.
 */
TreeMoments secondOrderStep(const TreeMoments &value, double step,
                            double size) {
  TreeMoments full = implicitStep(value, step, size);
  TreeMoments half =
      implicitStep(implicitStep(value, step / 2.0, size), step / 2.0, size);
  for (size_t k = 2; k < value.mean.size(); ++k) {
    half.mean[k] = 2.0 * half.mean[k] - full.mean[k];
    half.second[k] = 2.0 * half.second[k] - full.second[k];
  }
  return half;
}
} // namespace

/**
 * Calculate the mean and variance of total branch length without recombination.
 *
 * The R wrapper validates the inputs and converts population sizes and epoch
 * boundaries to coalescent units. The oldest population size applies forever,
 * so its exact constant-population moments initialize the calculation. Finite
 * epochs are then incorporated from oldest to youngest, ending at the present.
 * With no boundaries, the exact terminal moments are returned directly.
 *
 * Adaptive step doubling compares one second-order step with two half steps.
 * A step is accepted when the maximum absolute difference, divided by
 * tolerance*(1 + abs(fine)), is at most one across both moments and all states.
 * The two-half-step result is retained. Step-size adjustments use the cubic
 * local-error scaling of the second-order method and a safety factor of 0.8.
 * This is a local error criterion, not a bound on the final variance error.
 *
 * @param n Number of sampled haplotypes, at least two.
 * @param sizes Positive effective population sizes divided by N0, ordered
 *              from the present to the oldest epoch; the first entry is one.
 * @param times Strictly increasing positive boundary times divided by 2*N0.
 *              Must contain one fewer entry than sizes; may be empty.
 * @param tolerance Finite numerical tolerance strictly between zero and one.
 * @return Unnamed vector containing E[T] and Var(T), in coalescent units and
 *         squared coalescent units respectively. The R wrapper converts these
 *         to generations and maps them to segregating-site moments using theta.
 * @throws Rcpp exception if integration fails to converge, step sizes underflow,
 *         or calculated moments are invalid. User interrupts are checked between
 *         integration attempts. Input-domain validation belongs to the R wrapper.
 */
// [[Rcpp::export(name = ".segSitesTreeMoments", rng = false)]]
Rcpp::NumericVector segSitesTreeMoments(int n, Rcpp::NumericVector sizes,
                                        Rcpp::NumericVector times,
                                        double tolerance) {
  TreeMoments value(n);
  double size = sizes[sizes.size() - 1];
  double variance = 0.0;
  // The oldest epoch is constant forever, giving an exact terminal condition.
  for (int k = 2; k <= n; ++k) {
    // This is k times the mean waiting time: the branch-length contribution
    // of state k. Its variance equals its squared mean in the constant epoch.
    double waiting = 2.0 * size / (k - 1.0);
    value.mean[k] = value.mean[k - 1] + waiting;
    variance += waiting * waiting;
    value.second[k] = value.mean[k] * value.mean[k] + variance;
  }
  for (int epoch = times.size() - 1; epoch >= 0; --epoch) {
    double duration = times[epoch] - (epoch > 0 ? times[epoch - 1] : 0.0);
    double elapsed = 0.0;
    double step = duration / 10.0;
    int attempts = 0;
    while (elapsed < duration) {
      if (++attempts > 1000000)
        Rcpp::stop("Tree-moment integration did not converge");
      step = std::min(step, duration - elapsed);
      if (elapsed + step == elapsed)
        Rcpp::stop("Tree-moment integration step underflow");
      TreeMoments coarse = secondOrderStep(value, step, sizes[epoch]);
      TreeMoments fine =
          secondOrderStep(secondOrderStep(value, step / 2.0, sizes[epoch]),
                          step / 2.0, sizes[epoch]);
      double error = 0.0;
      for (int k = 2; k <= n; ++k) {
        error =
            std::max(error, std::abs(fine.mean[k] - coarse.mean[k]) /
                                (tolerance * (1.0 + std::abs(fine.mean[k]))));
        error =
            std::max(error, std::abs(fine.second[k] - coarse.second[k]) /
                                (tolerance * (1.0 + std::abs(fine.second[k]))));
      }
      if (!std::isfinite(error))
        Rcpp::stop("Non-finite tree moments during integration");
      if (error <= 1.0) {
        value = fine;
        elapsed += step;
      }
      step *= std::min(3.0, std::max(0.1, 0.8 * std::pow(std::max(error, 1e-12),
                                                         -1.0 / 3.0)));
      Rcpp::checkUserInterrupt();
    }
  }
  double mean = value.mean[n];
  variance = value.second[n] - mean * mean;
  if (!std::isfinite(mean) || !std::isfinite(variance) || variance < 0.0) {
    Rcpp::stop("Invalid calculated tree moments");
  }
  return Rcpp::NumericVector::create(mean, variance);
}
