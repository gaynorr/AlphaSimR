#ifndef MMECHOL_H
#define MMECHOL_H

// Helpers shared by the multiple trait solvers in mmeMultiTrait.cpp and
// mmeMultiKernel.cpp, which estimate unstructured covariance matrices among
// traits by Newton steps in their lower Cholesky factors.

#include "alphasimr.h"

// The parameters are the elements on and below the diagonal of Vu followed by
// those of Ve. The derivative of the covariance of a rotated row with respect
// to one of them is s_i*E (for Vu) or E (for Ve), where E has ones at (j,k)
// and (k,j).
inline void paramIndex(arma::uword t, arma::uvec& rowIdx, arma::uvec& colIdx){
  arma::uword nPar = t*(t+1)/2;
  rowIdx.set_size(nPar);
  colIdx.set_size(nPar);
  arma::uword p = 0;
  for(arma::uword k=0; k<t; ++k){
    for(arma::uword j=k; j<t; ++j){
      rowIdx(p) = j;
      colIdx(p) = k;
      ++p;
    }
  }
}

// Lower Cholesky factor of a covariance matrix that may be singular after an
// EM step on the boundary. A small ridge is added only when the plain
// factorization fails.
inline arma::mat lowerChol(const arma::mat& V){
  arma::mat L;
  if(arma::chol(L, V, "lower")){
    return L;
  }
  double ridge = 1.0e-10*std::max(arma::trace(V), 1.0e-300);
  arma::mat I = arma::eye(V.n_rows, V.n_cols);
  while(!arma::chol(L, V+ridge*I, "lower")){
    ridge *= 10.0;
  }
  return L;
}

// Jacobian of the elements of V = L*L' with respect to those of L, both
// listed by rowIdx and colIdx. With E_ab the unit matrix at (a,b),
// dV/dL_ab = E_ab*L' + L*E_ba, whose (j,k) element is
// [j==a]*L_kb + [k==a]*L_jb.
inline arma::mat cholJacobian(const arma::mat& L, const arma::uvec& rowIdx,
                       const arma::uvec& colIdx){
  arma::uword nPar = rowIdx.n_elem;
  arma::mat J(nPar, nPar, arma::fill::zeros);
  for(arma::uword q=0; q<nPar; ++q){
    arma::uword a = rowIdx(q), b = colIdx(q);
    for(arma::uword p=0; p<nPar; ++p){
      arma::uword j = rowIdx(p), k = colIdx(p);
      double val = 0.0;
      if(j==a) val += L(k,b);
      if(k==a) val += L(j,b);
      J(p,q) = val;
    }
  }
  return J;
}

// Second order term of the Hessian with respect to L. The second derivative
// of V = L*L' with respect to L_ab and L_cd is zero unless b equals d, when
// it is E_ac + E_ca. The term therefore pairs elements of L in the same
// column and takes the gradient with respect to the element (a,c) of V, twice
// over when a equals c. gradV holds that gradient on and below its diagonal.
inline arma::mat cholCurvature(const arma::mat& gradV, const arma::uvec& rowIdx,
                        const arma::uvec& colIdx){
  arma::uword nPar = rowIdx.n_elem;
  arma::mat K(nPar, nPar, arma::fill::zeros);
  for(arma::uword q=0; q<nPar; ++q){
    for(arma::uword r=0; r<nPar; ++r){
      if(colIdx(q)!=colIdx(r)) continue;
      arma::uword a = rowIdx(q), c = rowIdx(r);
      K(q,r) = (a==c) ? 2.0*gradV(a,a) : gradV(std::max(a,c), std::min(a,c));
    }
  }
  return K;
}

// The same gradient laid out on and below the diagonal of a matrix
inline arma::mat gradMatrix(const arma::vec& g, arma::uword offset, arma::uword t,
                     const arma::uvec& rowIdx, const arma::uvec& colIdx){
  arma::mat G(t, t, arma::fill::zeros);
  for(arma::uword p=0; p<rowIdx.n_elem; ++p){
    G(rowIdx(p), colIdx(p)) = g(offset+p);
  }
  return G;
}

// The responses standardized trait by trait. The Newton steps floor the
// eigenvalues of their Hessian with one threshold across every trait's
// coordinates, and a trait's curvature there scales with the inverse square
// of its units: a trait recorded in units 1e8 times larger would have its
// curvature floored and its updates damped, and the fit could stop before
// that trait had converged. Dividing each trait by the standard deviation of
// its residuals after the least squares fit puts every trait on the same
// scale, and mvResult maps the results back exactly.
struct TraitScale {
  arma::mat Y;     // the standardized responses
  arma::rowvec sd; // the divisor of each trait
};

inline TraitScale traitScale(const arma::mat& Ys){
  TraitScale ts;
  ts.sd = arma::stddev(Ys, 0, 0);
  for(arma::uword j=0; j<ts.sd.n_elem; ++j){
    if(!(ts.sd(j)>0.0) || !std::isfinite(ts.sd(j))){
      ts.sd(j) = 1.0;
    }
  }
  ts.Y = Ys;
  ts.Y.each_row() /= ts.sd;
  return ts;
}

#endif
