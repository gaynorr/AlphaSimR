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

// Columns of a fixed effects design X that are linearly independent, in
// increasing order.
//
// Every solver drops the rest rather than inverting a singular X'*V^-1*X.
// The dropped columns lie in the span of the kept ones, so the column space
// of X, and with it the REML likelihood, the variance components, the BLUPs
// and the fitted values, is unchanged. The dropped columns get a fixed effect
// of zero, which is one of the many solutions when X is not of full rank.
//
// Columns are taken greedily: a column is kept if its squared length after
// projecting out the columns already kept is more than a small fraction of
// its squared length. Which of a dependent set is dropped therefore depends
// on the order of consideration. The intercept, the first column whose
// entries are all the same nonzero value, is considered before every other
// column, so it is always kept and its estimate keeps its usual meaning
// whatever columns X also contains. The rest are considered in their order
// in X.
inline arma::uvec independentColumns(const arma::mat& X){
  arma::uword nCol = X.n_cols;
  std::vector<arma::uword> order;
  order.reserve(nCol);
  for(arma::uword j=0; j<nCol; ++j){
    if( (X.n_rows>0) && (X(0,j)!=0.0) && arma::all(X.col(j)==X(0,j)) ){
      order.push_back(j);
      break;
    }
  }
  for(arma::uword j=0; j<nCol; ++j){
    if(order.empty() || (order[0]!=j)){
      order.push_back(j);
    }
  }
  arma::mat XtX = X.t()*X;
  std::vector<arma::uword> keep;
  for(arma::uword j : order){
    double djj = XtX(j,j);
    if(!(djj>0.0)){
      continue;
    }
    double res = djj;
    if(!keep.empty()){
      arma::uvec kk(keep);
      arma::uvec jj(1);
      jj(0) = j;
      arma::vec c = XtX.submat(kk,jj);
      arma::vec b;
      // The columns kept so far are independent, so this system is
      // nonsingular and the solve only fails on rounding trouble, in which
      // case the column is treated as dependent
      if(!arma::solve(b, XtX.submat(kk,kk), c, arma::solve_opts::no_approx)){
        continue;
      }
      res = djj - arma::dot(c,b);
    }
    if(res > 1e-10*djj){
      keep.push_back(j);
    }
  }
  return arma::sort(arma::uvec(keep));
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

// Single trait RR-BLUP by FaST-LMM, defined in MME.cpp
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M);

#endif
