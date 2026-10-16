#ifndef MME_H
#define MME_H

// Helpers shared by the mixed model solvers in MME.cpp, mmeSingleKernel.cpp,
// mmeMultiTrait.cpp and mmeMultiKernel.cpp.

#include "alphasimr.h"

// Symmetric eigendecomposition through LAPACK's dsyevr, defined in MME.cpp.
// If calcVec is false, eigvec is not used.
int eigen2(arma::vec& eigval, arma::mat& eigvec, arma::mat X,
           bool calcVec = true);

// Chooses between the two ways of solving a single trait RR-BLUP model when
// the variance components are already known.
//
// Henderson's mixed model equations work with a square matrix of the fixed
// effects plus the markers. Forming M'M costs about n*m^2 and factoring the
// result about m^3/6, so the marker side costs roughly
//
//   n*m^2 + m^3/6
//
// The same model written on the records works with a square matrix of the
// records, because M*M' and M'M share their nonzero eigenvalues. Forming
// M*M' costs about n^2*m and factoring it about n^3/6, so the record side
// costs roughly
//
//   n^2*m + n^3/6
//
// One expression is the other with n and m exchanged, so the smaller of the
// two dimensions names the cheaper method and the crossover sits at n equal
// to m whatever the shared constant on the factorization turns out to be.
// Memory crosses over at the same place, since the two methods hold an n by n
// and an m by m matrix respectively on top of the genotypes they share. The
// costs are written out rather than reduced to a comparison of n against m so
// that the reasoning is visible, and so that there is somewhere to put a
// correction if the two sides ever stop being symmetric.
inline bool useRecordSide(arma::uword n, arma::uword m){
  const double dn = double(n), dm = double(m);
  const double recordSide = dn*dn*dm + dn*dn*dn/6.0;
  const double markerSide = dn*dm*dm + dm*dm*dm/6.0;
  return recordSide < markerSide;
}

// Inverts a symmetric positive definite matrix and returns its log
// determinant. One Cholesky factorization supplies both, so the determinant
// that a REML likelihood needs costs nothing beyond the inverse itself.
// Returns false if the matrix is not positive definite.
inline bool invAndLogDet(arma::mat& out, double& logDet, const arma::mat& in){
  arma::mat L;
  if(!arma::chol(L, in, "lower")){
    return false;
  }
  arma::mat Linv;
  if(!arma::inv(Linv, arma::trimatl(L))){
    return false;
  }
  logDet = 2.0*accu(log(L.diag()));
  out = Linv.t()*Linv;
  return true;
}

// Fit of the fixed effects at a given delta, shared by the FaST-LMM
// objective below and by the final solve.
struct FastLMMFit {
  arma::vec beta;      // generalized least squares solution
  double r = 0.0;      // (y-X*beta)'*H^-1*(y-X*beta)
  double logDetH = 0.0; // log|H|
  double logDetA = 0.0; // log|X'*H^-1*X|
  bool ok = false;
};

// The mixed model variance is V = Vu*(M*M'+delta*I) = Vu*H. Writing the
// spectral decomposition of M*M' as U*diag(s2)*U' together with a null space
// of dimension nNull, any quadratic form splits into a part inside the span
// of U and a part outside it,
//
//   a'*H^-1*b = sum_i at_i*bt_i/(s2_i+delta) + (a'*b-at'*bt)/delta
//
// for at = U'a and bt = U'b. The pieces outside the span do not depend on
// delta and are supplied as yy0, Xy0 and XX0, so one decomposition serves
// every value of delta and each evaluation is linear in the rank of M*M'.
inline FastLMMFit fastLMMFit(double delta, const arma::vec& yt,
                             const arma::mat& Xt, const arma::vec& s2,
                             double yy0, const arma::vec& Xy0,
                             const arma::mat& XX0, double nNull){
  FastLMMFit out;
  out.ok = false;
  arma::vec dInv = 1.0/(s2+delta);
  arma::mat XtD = Xt;
  XtD.each_col() %= dInv;
  arma::mat A = Xt.t()*XtD + XX0/delta;
  A = 0.5*(A+A.t()); // rounding only
  arma::vec b = XtD.t()*yt + Xy0/delta;
  double c = accu(yt%yt%dInv) + yy0/delta;
  arma::mat L;
  if(!arma::chol(L, A, "lower")){
    return out;
  }
  out.logDetA = 2.0*accu(log(L.diag()));
  out.beta = arma::solve(arma::trimatu(L.t()),
                         arma::solve(arma::trimatl(L), b));
  out.r = c - accu(b%out.beta);
  if(out.r<=0.0) out.r = arma::datum::eps;
  out.logDetH = accu(log(s2+delta)) + nNull*log(delta);
  out.ok = true;
  return out;
}

// FaST-LMM REML objective and its minimizer, defined in MME.cpp
Rcpp::List objREMLFastLMM(double param, Rcpp::List args);
double fastLMMDelta(Rcpp::List args);

// Single trait RR-BLUP by FaST-LMM, defined in MME.cpp
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M);

#endif
