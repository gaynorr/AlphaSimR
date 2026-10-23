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

// Chooses between the two ways of fitting a single trait REML model with
// several random terms by average information.
//
// The record side works with V, the n by n variance matrix of the records.
// Every likelihood evaluation forms, factors and inverts it and forms P, at
// about n^3, and the traces of the terms add about n^2*c.
//
// The coefficient side works with Henderson's mixed model equations, a c by
// c system for c fixed and random effects. Forming T'*T for the design T of
// all the effects costs n*c^2 once, after which every evaluation factors and
// inverts the coefficient matrix, at about c^3, and forms the average
// information from k+2 vectors of records, at about n*c*(k+2). extraCoef is
// any further one-off cost of putting the model in that form, such as
// decomposing the covariance matrices of the random terms.
//
// The one-off cost is set against a typical number of evaluations, counting
// step halvings. The crossover sits close to c equal to n, and memory
// crosses over at the same place, since the two sides hold n^2 and c^2
// doubles. The costs are written out so that the reasoning is visible.
inline bool useCoefficientSide(arma::uword n, arma::uword c, arma::uword k,
                               double extraCoef = 0.0){
  const double nEval = 10.0;
  const double dn = double(n), dc = double(c), dk = double(k);
  const double recordSide = nEval*(dn*dn*dn + dn*dn*dc);
  const double coefficientSide = dn*dc*dc + extraCoef +
    nEval*(dc*dc*dc + dn*dc*(dk+2.0));
  return coefficientSide < recordSide;
}

// Inverts a symmetric positive definite matrix and returns its log
// determinant. One Cholesky factorization supplies both, so the determinant
// that a REML likelihood needs costs nothing beyond the inverse itself.
// Returns false if the matrix is not positive definite.
// Eigenvalues of a cross product such as M'*M, M*M' or Z*K*Z' that are no
// larger than this are zero up to rounding. Every solver applies the same
// rule, to the eigenvalues of the cross product and relative to the number of
// records, whose sums carry the rounding. Whichever of M'*M and M*M' is
// decomposed, and whether or not M is padded with columns of zeros, the same
// directions are therefore held.
inline double rankCutoff(double maxEig, arma::uword nRecords){
  return double(nRecords)*arma::datum::eps*std::max(maxEig, 0.0);
}

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
  arma::mat L;
  if(!arma::chol(L, A, "lower")){
    return out;
  }
  out.logDetA = 2.0*accu(log(L.diag()));
  out.beta = arma::solve(arma::trimatu(L.t()),
                         arma::solve(arma::trimatl(L), b));
  // The residual quadratic form is summed from the residuals themselves
  // rather than found as y'*H^-1*y minus the part explained by beta, which
  // would subtract two large numbers and lose the residual variation when y
  // has a large mean. The part outside the span of U can only be reached
  // through the summaries, so callers keep it small by taking a least squares
  // fit off y first; see removeLeastSquares.
  arma::vec et = yt - Xt*out.beta;
  double r0 = (yy0 - 2.0*accu(Xy0%out.beta) +
               as_scalar(out.beta.t()*XX0*out.beta))/delta;
  if(r0<0.0) r0 = 0.0; // rounding only
  out.r = accu(et%et%dInv) + r0;
  if(out.r<=0.0) out.r = arma::datum::eps;
  out.logDetH = accu(log(s2+delta)) + nNull*log(delta);
  out.ok = true;
  return out;
}

// Columns of a fixed effects design X that are linearly independent, in
// increasing order.
//
// Every solver drops the rest rather than inverting a singular X'*V^-1*X.
// The dropped columns lie in the span of the kept ones, so the column space
// of X, and with it the REML likelihood, the variance components, the BLUPs
// and the fitted values, is unchanged. The dropped columns get a fixed effect
// of zero, which is one of the many solutions when X is not of full rank.
//
// Columns are taken greedily: a column is kept if its length after
// projecting out the columns already kept is more than a small fraction of
// its length. Which of a dependent set is dropped therefore depends on the
// order of consideration. The intercept, the first column whose entries are
// all the same nonzero value, is considered before every other column, so it
// is always kept and its estimate keeps its usual meaning whatever columns X
// also contains. The rest are considered in their order in X.
//
// The projection is done on the columns themselves, by Gram-Schmidt with a
// second pass to restore orthogonality, rather than through X'*X. Squaring
// the columns would lose a covariate whose variation is small next to its
// mean, such as 1e6+c(-1,0,1) beside an intercept, to rounding, although it
// is plainly independent; with the intercept projected out first, such a
// column is in effect centered before its length is judged. The cutoff
// allows for the rounding error of the projection, which grows with the
// number of records.
inline arma::uvec independentColumns(const arma::mat& X){
  const arma::uword n = X.n_rows;
  const arma::uword nCol = X.n_cols;
  std::vector<arma::uword> order;
  order.reserve(nCol);
  for(arma::uword j=0; j<nCol; ++j){
    if( (n>0) && (X(0,j)!=0.0) && arma::all(X.col(j)==X(0,j)) ){
      order.push_back(j);
      break;
    }
  }
  for(arma::uword j=0; j<nCol; ++j){
    if(order.empty() || (order[0]!=j)){
      order.push_back(j);
    }
  }
  const double tol = 100.0*double(std::max<arma::uword>(n, 1))*
    arma::datum::eps;
  std::vector<arma::vec> basis;
  std::vector<arma::uword> keep;
  for(arma::uword j : order){
    arma::vec v = X.col(j);
    double len = arma::norm(v);
    if(!(len>0.0) || !std::isfinite(len)){
      continue;
    }
    for(int pass=0; pass<2; ++pass){
      for(const arma::vec& q : basis){
        v -= arma::dot(q, v)*q;
      }
    }
    double res = arma::norm(v);
    if(res > tol*len){
      basis.push_back(v/res);
      keep.push_back(j);
    }
  }
  return arma::sort(arma::uvec(keep));
}

// The fixed effects design in the form the solvers fit.
//
// The independent columns of X are replaced by an orthonormal basis Q for
// their span, with X.cols(keep) = Q*R. REML depends on X only through its
// column space, so the variance components and random effects are the same
// with Q as with X, while Q is perfectly conditioned however X is scaled.
// Fitting with X directly would lose precision to a covariate with a large
// offset or scale, whose cross products dwarf the others. The fixed effects
// for X are R^-1 times those for Q, and log|X'*V^-1*X| is
// log|Q'*V^-1*Q| + 2*logDetR, so a likelihood found with Q is converted by
// subtracting logDetR for each trait.
struct FixedDesign {
  arma::uvec keep; // columns of X that are used
  arma::mat Q;     // orthonormal basis for their span
  arma::mat R;     // X.cols(keep) = Q*R
  double logDetR = 0.0;
};

inline FixedDesign fixedDesign(const arma::mat& X){
  // A column with a missing or infinite value has no finite length, and
  // would otherwise be dropped as if it depended on the others, leaving a
  // fit that looks successful without the covariate
  if(!X.is_finite()){
    Rcpp::stop("X must not contain missing or infinite values");
  }
  FixedDesign f;
  f.keep = independentColumns(X);
  if(f.keep.n_elem>0){
    if(!arma::qr_econ(f.Q, f.R, X.cols(f.keep))){
      Rcpp::stop("QR decomposition of X failed");
    }
    f.logDetR = arma::accu(arma::log(arma::abs(f.R.diag())));
  }
  return f;
}

// Fixed effects for the columns of X from those estimated for Q, one column
// per trait, with zeros for the columns that were dropped
inline arma::mat fixedEffectsFromBasis(const FixedDesign& f,
                                       const arma::mat& Bq,
                                       arma::uword nCol){
  arma::mat out(nCol, Bq.n_cols, arma::fill::zeros);
  out.rows(f.keep) = arma::solve(arma::trimatu(f.R), Bq);
  return out;
}

// Takes a least squares fit on the fixed effects off the responses.
//
// REML, and the solution for everything but the fixed effects, are unchanged
// when Y is replaced by Y-X*B0 for any B0, while the fixed effects move by
// exactly -B0. With B0 the least squares fit, the responses are left on the
// scale of their own variation, so a large mean cannot swamp the sums of
// squares the likelihood is built from. Callers add B0 back to the fixed
// effects they estimate. X must have full column rank.
inline arma::mat removeLeastSquares(const arma::mat& X, const arma::mat& Y,
                                    arma::mat& B0){
  if(!arma::solve(B0, X, Y)){
    B0.zeros(X.n_cols, Y.n_cols);
    return Y;
  }
  return Y - X*B0;
}

// Whether a symmetric matrix is positive definite, by attempting its Cholesky
// factorization. That costs about a sixth of an eigendecomposition, and the
// covariance matrices passed to the solvers, such as relationship matrices
// with a small ridge, are usually positive definite, so psdKernel tries
// it first and only decompose a matrix it fails on.
inline bool choleskySucceeds(const arma::mat& Ks){
  arma::mat L;
  return arma::chol(L, Ks, "lower");
}

// K as the positive semidefinite covariance matrix that is both fitted and
// used for the BLUPs. Stops with msg if K has a clearly negative eigenvalue.
// Otherwise K is returned symmetrized, except that negative eigenvalues at
// rounding level, which the likelihood treats as zero, are set to zero in
// it too. Without that the BLUPs, which are computed with K itself, would
// describe a slightly different model from the one fitted, and a small
// residual variance would amplify the difference. A positive definite K,
// recognized by its Cholesky factorization, has no such eigenvalue and is
// returned without any decomposition. Otherwise the eigenvalues are
// computed, and the full decomposition only when such an eigenvalue is
// present.
inline arma::mat psdKernel(const arma::mat& K, const std::string& msg){
  arma::mat Ks = 0.5*(K+K.t());
  if(Ks.n_rows==0){
    return Ks;
  }
  if(choleskySucceeds(Ks)){
    return Ks;
  }
  arma::vec lambda;
  if(!arma::eig_sym(lambda, Ks)){
    Rcpp::stop("eigendecomposition of K failed");
  }
  double maxAbs = arma::abs(lambda).max();
  if(lambda.min() < -std::sqrt(arma::datum::eps)*maxAbs){
    Rcpp::stop(msg);
  }
  if(lambda.min()<0.0){
    arma::mat Q;
    if(!arma::eig_sym(lambda, Q, Ks)){
      Rcpp::stop("eigendecomposition of K failed");
    }
    lambda.elem(find(lambda<0.0)).zeros();
    Ks = Q*arma::diagmat(lambda)*Q.t();
    Ks = 0.5*(Ks+Ks.t()); // rounding only
  }
  return Ks;
}

// Z*K*Z', the covariance on the records implied by a covariance K among m
// levels and an n by m design Z.
//
// When Z is an incidence matrix, with a single one in every row and zeros
// elsewhere, as it is whenever each record belongs to exactly one level,
// Z*K*Z' is K with its rows and columns picked out by the level of each
// record. That costs n^2 against about n*m^2 for the dense products, which
// matters most when K covers many more levels than there are records, such
// as individuals with no record that are to be predicted. Checking Z costs
// n*m, no more than reading it. Any other Z falls back to the products.
inline arma::mat kernelOnRecords(const arma::mat& Z, const arma::mat& K){
  const arma::uword n = Z.n_rows;
  const arma::uword m = Z.n_cols;
  arma::uvec lev(n, arma::fill::zeros);
  std::vector<bool> seen(n, false);
  bool incidence = true;
  for(arma::uword j=0; (j<m) && incidence; ++j){
    const double* col = Z.colptr(j);
    for(arma::uword i=0; i<n; ++i){
      if(col[i]==0.0){
        continue;
      }
      if( (col[i]!=1.0) || seen[i] ){
        incidence = false;
        break;
      }
      seen[i] = true;
      lev(i) = j;
    }
  }
  if(incidence){
    for(arma::uword i=0; i<n; ++i){
      if(!seen[i]){
        incidence = false;
        break;
      }
    }
  }
  if(incidence){
    return K.submat(lev, lev);
  }
  return Z*K*Z.t();
}

// Places the rows of B, estimated on the columns keep of a fixed effects
// design, into a matrix with one row per original column. Rows for dropped
// columns are zero.
inline arma::mat expandFixedEffects(const arma::mat& B, const arma::uvec& keep,
                                    arma::uword nCol){
  arma::mat out(nCol, B.n_cols, arma::fill::zeros);
  out.rows(keep) = B;
  return out;
}

// FaST-LMM REML objective and its minimizer, defined in MME.cpp
Rcpp::List objREMLFastLMM(double param, Rcpp::List args);
double fastLMMDelta(Rcpp::List args);

// Multiple trait REML on a kernel among the records, defined in
// mmeMultiTrait.cpp
void remlMultiTrait(const arma::mat& Y, const arma::mat& X,
                    const arma::mat& F, const arma::mat& Kn, bool useFactor,
                    int maxIter, double tol, arma::mat& Vu, arma::mat& Ve);

// Pieces of the RR-BLUP models shared with the fast mode of the multikernel
// models in mmeFastMK.cpp, defined in MME.cpp
arma::mat makeX(arma::uvec& x);
arma::Mat<unsigned char> packChipPlanes(const PackedGeno& g, const LociPos& pos,
                                        int nThreads);
arma::rowvec negMeanProduct(const arma::rowvec& Mmean, const arma::mat& U);
void domAverageEffects(const arma::mat& freqMat, const arma::uvec& fixed,
                       const arma::mat& ua, const arma::mat& ud,
                       const arma::rowvec& meanDAll, arma::uword ploidy,
                       int nThreads, arma::mat& alpha, arma::mat& d);
void scaAverageEffects(const arma::rowvec& Mmean1, const arma::rowvec& Mmean2,
                       const arma::rowvec& p12, const arma::uvec& fixed,
                       const arma::mat& u0, const arma::mat& u1,
                       const arma::mat& u2, const arma::rowvec& meanDAll,
                       int nThreads, arma::mat& alpha1, arma::mat& alpha2,
                       arma::mat& d);

// Single trait RR-BLUP by FaST-LMM, defined in MME.cpp
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M);

#endif
