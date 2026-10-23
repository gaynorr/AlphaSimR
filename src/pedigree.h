#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <RcppArmadillo.h>
#include <vector>

// A pedigree in the form the relationship functions work on: rows in an
// order where every parent comes before its progeny, 0-based, with -1 for
// an unknown parent. A doubled haploid has its one parent as both mother
// and father. An inbred founder is marked as a doubled haploid with no
// parents. All individuals have the same ploidy, an even number.
struct Ped{
  arma::uword n;
  std::vector<int> mother;
  std::vector<int> father;
  std::vector<bool> isDH;
  int ploidy;
};

// Builds a Ped from the 1-based vectors that R passes, with 0 for an
// unknown parent, and checks that it is in a usable order.
Ped readPed(const Rcpp::IntegerVector& mother,
            const Rcpp::IntegerVector& father,
            const Rcpp::IntegerVector& isDH,
            int ploidy);

// Inbreeding coefficients (F) and Mendelian sampling variances (D) of every
// individual in the pedigree. knownF is empty, or holds the F already
// known with NaN for the rest, as from the SimParam's cache.
void pedInbreeding(const Ped& ped, const arma::vec& knownF, arma::vec& F,
                   arma::vec& D);

// The inverse of the numerator relationship matrix, written so that it
// exists even when the matrix itself is singular.
//
// An individual whose Mendelian sampling variance is zero, such as a
// doubled haploid of a fully inbred parent or the F1 of two fully inbred
// parents, has a breeding value that is an exact linear function of its
// parents', and the relationship matrix is then singular. Such individuals
// are left out of the basis and expressed through W, so that the breeding
// values of the whole pedigree are u = W*b, where b ~ N(0, Ainv^-1) covers
// the basis individuals only. W has a single 1 in each basis individual's
// row.
struct PedAinv{
  arma::sp_mat Ainv; // nBasis by nBasis
  arma::sp_mat W;    // n by nBasis
  arma::uvec basis;  // 0-based pedigree rows of the basis individuals
};

PedAinv buildPedAinv(const Ped& ped, const arma::vec& D);

#endif
