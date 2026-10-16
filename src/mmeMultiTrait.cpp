#include "mme.h"

// Multiple trait, single kernel mixed models,
//
//   Y = X*B + G + E,  vec(G) ~ N(0, Vu (x) Kn),  vec(E) ~ N(0, Ve (x) I),
//
// for t traits on n records, with Kn the n by n kernel on the records and Vu
// and Ve unstructured t by t covariance matrices. The variance components are
// estimated by REML. This is an independent implementation from the published
// description of Zhou and Stephens (2014), doi:10.1038/nmeth.2848, together
// with textbook REML algebra for the EM and average information steps; no code
// is taken from any other software.
//
// The method rests on two rotations. Writing Kn = U*diag(s)*U', the rotated
// records U'*Y are independent across rows, row i having covariance
// s_i*Vu + Ve. Writing Ve = L*L' and L^-1*Vu*L^-T = Q*diag(lambda)*Q', the
// matrix T = Q'*L^-1 then diagonalizes every row at once, because
//
//   T*(s_i*Vu + Ve)*T' = diag(1 + s_i*lambda).
//
// B is unrestricted, so in the coordinates given by T the model separates into
// t univariate weighted least squares problems that share the rotation, and
// every quantity the likelihood and its derivatives need is linear in n.

namespace {

// The rotated problem. Only the rows of U'*Y with a positive eigenvalue need to
// be held when the kernel is of low rank, because every row of the null space
// has the same covariance, Ve. Those rows enter only through the sums of
// squares and products below, which are zero when every row is held.
struct MvData {
  arma::mat Yt;   // rotated responses, n1 by t
  arma::mat Xt;   // rotated fixed effects design, n1 by c
  arma::vec s;    // eigenvalues of the kernel for the held rows, n1
  arma::mat YY0;  // Y0'*Y0 over the rows that are not held, t by t
  arma::mat XY0;  // X0'*Y0, c by t
  arma::mat XX0;  // X0'*X0, c by c
  double nNull = 0.0; // number of rows that are not held
  double n = 0.0;     // number of records
  double r = 0.0;     // number of positive eigenvalues
};

// Everything known at one value of (Vu, Ve). The quantities in transformed
// coordinates are kept so that the average information matrix can be formed
// from them without repeating the fit.
struct MvEval {
  bool ok = false;
  double LL = 0.0;
  arma::mat B;      // fixed effects, c by t
  arma::mat Su;     // sum_i s_i*(a_i*a_i' - P_ii), t by t
  arma::mat Se;     // sum_i (a_i*a_i' - P_ii), t by t
  arma::mat T;      // Q'*L^-1
  arma::mat D;      // 1/(1+s_i*lambda_k), n1 by t
  arma::mat Ahat;   // P*y for the held rows, transformed coordinates
  arma::mat Q0;     // sum of Ahat*Ahat' over the null rows
  arma::mat XA0;    // X0'*Ahat over the null rows, c by t
  arma::field<arma::mat> Winv; // inverse of X'*V^-1*X for each transformed trait
};

// Restricted log-likelihood, fixed effects and scores at (Vu, Ve).
//
// In transformed coordinates the generalized least squares problem for trait k
// has weights d_ik = 1/(1+s_i*lambda_k) on the held rows and 1 on the null
// rows. With a_i the rows of P*y and P_ii the diagonal blocks of P, the
// derivative of the restricted log-likelihood with respect to Vu is
// 0.5*sum_i s_i*(a_i*a_i' - P_ii), and with respect to Ve the same sum without
// s_i. Both are computed in transformed coordinates, where P_ii is diagonal,
// and returned in the original ones.
MvEval evalMv(const MvData& d, const arma::mat& Vu, const arma::mat& Ve){
  MvEval out;
  const arma::uword t = d.Yt.n_cols;
  const arma::uword c = d.Xt.n_cols;
  arma::mat L;
  if(!arma::chol(L, Ve, "lower")){
    return out;
  }
  arma::mat Linv;
  if(!arma::inv(Linv, arma::trimatl(L))){
    return out;
  }
  arma::mat A = Linv*Vu*Linv.t();
  A = 0.5*(A+A.t());
  arma::vec lambda;
  arma::mat Q;
  if(!arma::eig_sym(lambda, Q, A)){
    return out;
  }
  // Vu may sit on the boundary of the parameter space, where a rounding error
  // leaves a slightly negative eigenvalue. A clearly negative one means Vu is
  // not a covariance matrix.
  double lamMax = std::max(lambda.max(), 1.0);
  if(lambda.min() < -1.0e-8*lamMax){
    return out;
  }
  lambda.elem(find(lambda<0.0)).zeros();
  out.T = Q.t()*Linv;
  const arma::mat& T = out.T;

  arma::mat Yh = d.Yt*T.t();
  arma::mat XYh0 = d.XY0*T.t();
  arma::mat YYh0 = T*d.YY0*T.t();
  out.D = 1.0/(1.0+d.s*lambda.t());
  const arma::mat& D = out.D;

  arma::mat Bh(c, t);
  arma::mat H(d.Xt.n_rows, t); // leverages x_i'*W_k^-1*x_i
  arma::vec trWX0(t);
  out.Winv.set_size(t);
  double logDetW = 0.0;
  for(arma::uword k=0; k<t; ++k){
    arma::mat XtD = d.Xt;
    XtD.each_col() %= D.col(k);
    arma::mat W = d.Xt.t()*XtD + d.XX0;
    W = 0.5*(W+W.t());
    double ld;
    if(!invAndLogDet(out.Winv(k), ld, W)){
      Rcpp::stop("X is not of full column rank");
    }
    logDetW += ld;
    Bh.col(k) = out.Winv(k)*(XtD.t()*Yh.col(k) + XYh0.col(k));
    H.col(k) = sum((d.Xt*out.Winv(k))%d.Xt, 1);
    trWX0(k) = accu(out.Winv(k)%d.XX0);
  }

  arma::mat R = Yh - d.Xt*Bh;
  out.Ahat = D%R;
  out.Q0 = YYh0 - XYh0.t()*Bh - Bh.t()*XYh0 + Bh.t()*d.XX0*Bh;
  out.Q0 = 0.5*(out.Q0+out.Q0.t());
  out.XA0 = XYh0 - d.XX0*Bh;
  arma::mat Pd = D - D%D%H; // diagonals of the blocks P_ii

  double yPy = accu(out.Ahat%R) + trace(out.Q0);
  double logDetVe = 2.0*accu(log(L.diag()));
  double logDetV = d.n*logDetVe + accu(log(1.0+d.s*lambda.t()));
  double logDetXVX = logDetW - double(c)*logDetVe;
  double df = d.n*double(t) - double(c*t);
  out.LL = -0.5*(df*log(2.0*arma::datum::pi) + logDetV + logDetXVX + yPy);

  arma::mat AhatS = out.Ahat;
  AhatS.each_col() %= d.s;
  arma::mat Suh = out.Ahat.t()*AhatS - arma::diagmat(d.s.t()*Pd);
  arma::mat Seh = out.Ahat.t()*out.Ahat - arma::diagmat(sum(Pd, 0)) +
    out.Q0 - arma::diagmat(d.nNull - trWX0);
  out.Su = T.t()*Suh*T;
  out.Su = 0.5*(out.Su+out.Su.t());
  out.Se = T.t()*Seh*T;
  out.Se = 0.5*(out.Se+out.Se.t());

  // Back to the original coordinates, B = Bh*(T')^-1 = Bh*Q'*L'
  out.B = Bh*Q.t()*L.t();
  out.ok = true;
  return out;
}

// The parameters are the elements on and below the diagonal of Vu followed by
// those of Ve. The derivative of the covariance of a rotated row with respect
// to one of them is s_i*E (for Vu) or E (for Ve), where E has ones at (j,k)
// and (k,j).
void paramIndex(arma::uword t, arma::uvec& rowIdx, arma::uvec& colIdx){
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

// The refinement steps work with the lower Cholesky factors of Vu and Ve
// rather than with the matrices themselves. Any value of the factors gives a
// valid pair of covariance matrices, so a step cannot leave the parameter
// space, and an estimate whose correlation approaches one, which is common
// with few records, is approached smoothly instead of being cut off.
arma::vec packChol(const arma::mat& Lu, const arma::mat& Le,
                   const arma::uvec& rowIdx, const arma::uvec& colIdx){
  arma::uword nPar = rowIdx.n_elem;
  arma::vec phi(2*nPar);
  for(arma::uword p=0; p<nPar; ++p){
    phi(p) = Lu(rowIdx(p), colIdx(p));
    phi(nPar+p) = Le(rowIdx(p), colIdx(p));
  }
  return phi;
}

void unpackChol(const arma::vec& phi, arma::mat& Vu, arma::mat& Ve,
                arma::uword t, const arma::uvec& rowIdx,
                const arma::uvec& colIdx){
  arma::uword nPar = rowIdx.n_elem;
  arma::mat Lu(t, t, arma::fill::zeros), Le(t, t, arma::fill::zeros);
  for(arma::uword p=0; p<nPar; ++p){
    Lu(rowIdx(p), colIdx(p)) = phi(p);
    Le(rowIdx(p), colIdx(p)) = phi(nPar+p);
  }
  Vu = Lu*Lu.t();
  Ve = Le*Le.t();
}

// Lower Cholesky factor of a covariance matrix that may be singular after an
// EM step on the boundary. A small ridge is added only when the plain
// factorization fails.
arma::mat lowerChol(const arma::mat& V){
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
arma::mat cholJacobian(const arma::mat& L, const arma::uvec& rowIdx,
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
arma::mat cholCurvature(const arma::mat& gradV, const arma::uvec& rowIdx,
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

// Gradient of the restricted log-likelihood, 0.5*trace(S*E) for each
// parameter.
arma::vec gradParam(const MvEval& ev, const arma::uvec& rowIdx,
                    const arma::uvec& colIdx){
  arma::uword nPar = rowIdx.n_elem;
  arma::vec g(2*nPar);
  for(arma::uword p=0; p<nPar; ++p){
    double w = (rowIdx(p)==colIdx(p)) ? 0.5 : 1.0;
    g(p) = w*ev.Su(rowIdx(p), colIdx(p));
    g(nPar+p) = w*ev.Se(rowIdx(p), colIdx(p));
  }
  return g;
}

// The same gradient laid out on and below the diagonal of a matrix
arma::mat gradMatrix(const arma::vec& g, arma::uword offset, arma::uword t,
                     const arma::uvec& rowIdx, const arma::uvec& colIdx){
  arma::mat G(t, t, arma::fill::zeros);
  for(arma::uword p=0; p<rowIdx.n_elem; ++p){
    G(rowIdx(p), colIdx(p)) = g(offset+p);
  }
  return G;
}

// Average information matrix of Gilmour et al. (1995), with elements
// 0.5*w_a'*P*w_b for the working variates w_a = dV/dtheta_a*P*y.
//
// In transformed coordinates the working variate of a parameter is
// s_i^e*M*a_i on held row i, with M = T*E*T' and e equal to 1 for Vu and 0 for
// Ve, and M*a_i on a null row for Ve. P is the diagonal weighting by d_ik less
// a projection that is separate for each transformed trait, so both terms of
// w_a'*P*w_b reduce to the weighted cross products of Ahat collected below,
// and no vector of length n*t is formed for any parameter.
arma::mat aiMatrix(const MvData& d, const MvEval& ev,
                   const arma::uvec& rowIdx, const arma::uvec& colIdx){
  const arma::uword t = d.Yt.n_cols;
  const arma::uword nPar = rowIdx.n_elem;
  const arma::uword nTot = 2*nPar;

  // C(e)_k = Ahat'*diag(d_k%s^e)*Ahat and G(e)_k = X'*diag(d_k%s^e)*Ahat
  arma::field<arma::mat> C(3, t), G(2, t);
  for(arma::uword k=0; k<t; ++k){
    arma::vec w = ev.D.col(k);
    for(arma::uword e=0; e<3; ++e){
      arma::mat Aw = ev.Ahat;
      Aw.each_col() %= w;
      C(e,k) = ev.Ahat.t()*Aw;
      if(e<2){
        G(e,k) = d.Xt.t()*Aw;
      }
      w %= d.s;
    }
  }

  // M for each parameter, and its exponent of s
  arma::field<arma::mat> M(nTot);
  arma::uvec ex(nTot);
  for(arma::uword p=0; p<nPar; ++p){
    arma::mat E(t, t, arma::fill::zeros);
    E(rowIdx(p), colIdx(p)) = 1.0;
    E(colIdx(p), rowIdx(p)) = 1.0;
    M(p) = ev.T*E*ev.T.t();
    M(nPar+p) = M(p);
    ex(p) = 1;
    ex(nPar+p) = 0;
  }

  // g_a for each transformed trait k, which is X'*V^-1*w_a for that trait
  arma::field<arma::mat> g(nTot);
  for(arma::uword a=0; a<nTot; ++a){
    g(a).set_size(d.Xt.n_cols, t);
    for(arma::uword k=0; k<t; ++k){
      g(a).col(k) = G(ex(a),k)*M(a).col(k);
    }
    if(ex(a)==0){
      g(a) += ev.XA0*M(a);
    }
  }

  arma::mat AI(nTot, nTot);
  for(arma::uword a=0; a<nTot; ++a){
    for(arma::uword b=0; b<=a; ++b){
      double val = 0.0;
      for(arma::uword k=0; k<t; ++k){
        val += as_scalar(M(a).row(k)*C(ex(a)+ex(b),k)*M(b).col(k));
        val -= as_scalar(g(a).col(k).t()*ev.Winv(k)*g(b).col(k));
      }
      if( (ex(a)==0) && (ex(b)==0) ){
        val += accu((M(a)*M(b))%ev.Q0.t());
      }
      AI(a,b) = AI(b,a) = 0.5*val;
    }
  }
  return AI;
}

bool isPosDef(const arma::mat& V){
  arma::mat L;
  return arma::chol(L, V, "lower");
}

// Rotates the problem by the eigendecomposition of the kernel F*F'.
//
// When F has fewer columns than rows, its thin singular value decomposition
// gives the eigenvectors with a positive eigenvalue directly and no n by n
// matrix is formed. The rest of the records then lie in the null space and
// are summarized. Otherwise Kn is decomposed itself, every row is held and an
// eigenvalue that is zero up to rounding is set to zero.
//
// Ur receives the eigenvectors of the held rows, n by n1.
MvData rotateData(const arma::mat& Y, const arma::mat& X, const arma::mat& F,
                  const arma::mat& Kn, bool useFactor, arma::mat& Ur){
  MvData d;
  const arma::uword n = Y.n_rows;
  d.n = double(n);
  if(useFactor){
    arma::mat Ufull, V;
    arma::vec sv;
    if(!arma::svd_econ(Ufull, sv, V, F, "left")){
      Rcpp::stop("Decomposition of the kernel failed");
    }
    double svMax = sv.n_elem>0 ? sv.max() : 0.0;
    double thr = svMax*double(std::max(F.n_rows, F.n_cols))*arma::datum::eps;
    arma::uvec keep = find(sv>thr);
    Ur = Ufull.cols(keep);
    Ufull.reset();
    d.s = square(sv(keep));
    d.Yt = Ur.t()*Y;
    d.Xt = Ur.t()*X;
    d.nNull = double(n) - double(keep.n_elem);
    // Subtracting the held part leaves the projection on the null space
    arma::mat Y0 = Y - Ur*d.Yt;
    arma::mat X0 = X - Ur*d.Xt;
    d.YY0 = Y0.t()*Y0;
    d.XY0 = X0.t()*Y0;
    d.XX0 = X0.t()*X0;
  }else{
    arma::vec eigval(n);
    Ur.set_size(n, n);
    if(eigen2(eigval, Ur, 0.5*(Kn+Kn.t()))!=0){
      Rcpp::stop("Decomposition of the kernel failed");
    }
    double eigMax = std::max(eigval.max(), 0.0);
    eigval.elem(find(eigval<=(double(n)*arma::datum::eps*eigMax))).zeros();
    d.s = eigval;
    d.Yt = Ur.t()*Y;
    d.Xt = Ur.t()*X;
    d.YY0.zeros(Y.n_cols, Y.n_cols);
    d.XY0.zeros(X.n_cols, Y.n_cols);
    d.XX0.zeros(X.n_cols, X.n_cols);
  }
  d.r = double(accu(d.s>0.0));
  if(d.r==0.0){
    Rcpp::stop("The kernel has no positive eigenvalue");
  }
  return d;
}

// Starting values. Each trait is first fitted on its own by the FaST-LMM REML
// search of Lippert et al. (2011), which supplies the diagonals of Vu and Ve.
// The off-diagonals take the correlation among the least squares residuals of
// the traits, which is the phenotypic correlation net of the fixed effects,
// for both matrices. If that leaves either matrix singular the off-diagonals
// start at zero instead.
void startValues(const MvData& d, arma::mat& Vu, arma::mat& Ve){
  const arma::uword t = d.Yt.n_cols;
  const arma::uword c = d.Xt.n_cols;
  arma::vec vu(t), ve(t);
  for(arma::uword k=0; k<t; ++k){
    arma::vec yt = d.Yt.col(k);
    arma::vec Xy0 = d.XY0.col(k);
    Rcpp::List args = Rcpp::List::create(Rcpp::Named("df")=d.n-double(c),
                                         Rcpp::Named("yt")=yt,
                                         Rcpp::Named("Xt")=d.Xt,
                                         Rcpp::Named("s2")=d.s,
                                         Rcpp::Named("yy0")=d.YY0(k,k),
                                         Rcpp::Named("Xy0")=Xy0,
                                         Rcpp::Named("XX0")=d.XX0,
                                         Rcpp::Named("nNull")=d.nNull);
    double delta = fastLMMDelta(args);
    FastLMMFit fit = fastLMMFit(delta, yt, d.Xt, d.s, d.YY0(k,k), Xy0,
                                d.XX0, d.nNull);
    if(!fit.ok){
      Rcpp::stop("X is not of full column rank");
    }
    vu(k) = fit.r/(d.n-double(c));
    ve(k) = delta*vu(k);
  }
  // A component estimated at the edge of the search would hold EM near zero
  // for many iterations, so each starts at no less than a thousandth of the
  // total. The kernel's scale enters through its mean eigenvalue.
  double sMean = accu(d.s)/d.n;
  for(arma::uword k=0; k<t; ++k){
    double tot = vu(k)*sMean + ve(k);
    if(!(tot>0.0)) tot = 1.0;
    vu(k) = std::max(vu(k), 1.0e-3*tot/std::max(sMean, arma::datum::eps));
    ve(k) = std::max(ve(k), 1.0e-3*tot);
  }

  // Correlation of the least squares residuals, from the rotated data
  arma::mat XX = d.Xt.t()*d.Xt + d.XX0;
  arma::mat XY = d.Xt.t()*d.Yt + d.XY0;
  arma::mat YY = d.Yt.t()*d.Yt + d.YY0;
  arma::mat Bols;
  arma::mat Rr;
  if(arma::solve(Bols, XX, XY, arma::solve_opts::no_approx)){
    Rr = YY - XY.t()*Bols;
  }else{
    Rr = arma::eye(t, t);
  }
  arma::vec sd = sqrt(arma::clamp(Rr.diag(), arma::datum::eps,
                                  arma::datum::inf));
  arma::mat rho = Rr/(sd*sd.t());
  rho.diag().ones();

  Vu = arma::diagmat(sqrt(vu))*rho*arma::diagmat(sqrt(vu));
  Ve = arma::diagmat(sqrt(ve))*rho*arma::diagmat(sqrt(ve));
  if(!isPosDef(Vu) || !isPosDef(Ve)){
    Vu = arma::diagmat(vu);
    Ve = arma::diagmat(ve);
  }
}

// REML estimation of Vu and Ve.
//
// The first few iterations are EM steps, which increase the likelihood from
// any start but converge slowly, as Zhou and Stephens (2014) describe. Newton
// steps follow, taken in the Cholesky factors of Vu and Ve. Their Hessian is
// the average information matrix carried through the Jacobian of V = L*L',
// plus the exact second order term of that map, which keeps the curvature
// where a factor approaches singularity. A direction of negative curvature
// has the sign of its eigenvalue reversed, so the step is always one of
// ascent. A step is halved up to ten times until the likelihood does not
// fall, and an EM step is taken instead if that fails.
//
// The EM step is the REML expectation-maximization update for a random term
// with r levels and a residual with n,
//
//   Vu <- Vu + Vu*Su*Vu/r,  Ve <- Ve + Ve*Se*Ve/n,
//
// written with the score matrices above. The expected sum of squares of the
// random term divides by s_i on each row, and this form cancels it, so rows
// with a small eigenvalue cause no loss of precision and rows with a zero one
// drop out.
//
// Iteration stops when the restricted log-likelihood changes by less than tol
// in a step taken after the initial EM steps.
MvEval fitMv(const MvData& d, arma::mat& Vu, arma::mat& Ve, int maxIter,
             double tol, int& iter){
  const arma::uword t = d.Yt.n_cols;
  const int nEM = 5;
  arma::uvec rowIdx, colIdx;
  paramIndex(t, rowIdx, colIdx);
  const arma::uword nPar = rowIdx.n_elem;
  const arma::uword nTot = 2*nPar;

  startValues(d, Vu, Ve);
  MvEval cur = evalMv(d, Vu, Ve);
  if(!cur.ok){
    Rcpp::stop("Starting values are not valid covariance matrices");
  }
  bool converged = false;
  for(iter=1; iter<=maxIter; ++iter){
    arma::mat VuNew, VeNew;
    MvEval next;
    bool accepted = false;
    if(iter>nEM){
      arma::mat Lu = lowerChol(Vu), Le = lowerChol(Ve);
      arma::mat AI = aiMatrix(d, cur, rowIdx, colIdx);
      arma::vec gV = gradParam(cur, rowIdx, colIdx);
      arma::mat J(nTot, nTot, arma::fill::zeros);
      J.submat(0, 0, nPar-1, nPar-1) = cholJacobian(Lu, rowIdx, colIdx);
      J.submat(nPar, nPar, nTot-1, nTot-1) = cholJacobian(Le, rowIdx, colIdx);
      arma::mat Kc(nTot, nTot, arma::fill::zeros);
      Kc.submat(0, 0, nPar-1, nPar-1) =
        cholCurvature(gradMatrix(gV, 0, t, rowIdx, colIdx), rowIdx, colIdx);
      Kc.submat(nPar, nPar, nTot-1, nTot-1) =
        cholCurvature(gradMatrix(gV, nPar, t, rowIdx, colIdx), rowIdx, colIdx);
      arma::vec gL = J.t()*gV;
      // Negative Hessian of the likelihood with respect to the factors
      arma::mat negH = J.t()*AI*J - Kc;
      negH = 0.5*(negH+negH.t());
      arma::vec hVal;
      arma::mat hVec;
      if(arma::eig_sym(hVal, hVec, negH)){
        hVal = arma::abs(hVal);
        double floorVal = 1.0e-10*std::max(hVal.max(), 1.0e-300);
        hVal.elem(find(hVal<floorVal)).fill(floorVal);
        arma::vec step = hVec*((hVec.t()*gL)/hVal);
        arma::vec phi = packChol(Lu, Le, rowIdx, colIdx);
        double alpha = 1.0;
        for(int h=0; h<10; ++h){
          unpackChol(phi+alpha*step, VuNew, VeNew, t, rowIdx, colIdx);
          next = evalMv(d, VuNew, VeNew);
          if(next.ok && (next.LL>=cur.LL)){
            accepted = true;
            break;
          }
          alpha *= 0.5;
        }
      }
    }
    if(!accepted){
      VuNew = Vu + Vu*cur.Su*Vu/d.r;
      VuNew = 0.5*(VuNew+VuNew.t());
      VeNew = Ve + Ve*cur.Se*Ve/d.n;
      VeNew = 0.5*(VeNew+VeNew.t());
      next = evalMv(d, VuNew, VeNew);
      if(!next.ok){
        Rcpp::stop("REML iteration left the parameter space");
      }
    }
    double change = next.LL - cur.LL;
    Vu = VuNew;
    Ve = VeNew;
    cur = next;
    if( (iter>nEM) && (std::abs(change)<tol) ){
      converged = true;
      break;
    }
  }
  if(!converged){
    iter = maxIter;
    Rf_warning("REML did not converge within maxIter iterations");
  }
  return cur;
}

// P*y on the records, n by t, in the original trait coordinates, for the
// residuals r = y - X*B at the final estimates.
//
// The fit works through T, which is ill conditioned when Ve is close to
// singular. That happens when the REML estimate of a residual correlation is
// at or near one, and the BLUPs would lose precision there. Each held row is
// therefore solved against its own covariance s_i*Vu + Ve, which stays well
// conditioned for a positive s_i whenever Vu is positive definite. A held
// row with a zero eigenvalue, and every null row, has covariance Ve itself,
// and the null rows are the projection of the residuals on the null space.
arma::mat recordPy(const arma::mat& Y, const arma::mat& X, const arma::mat& Ur,
                   const MvData& d, const MvEval& ev, const arma::mat& Vu,
                   const arma::mat& Ve){
  arma::mat R = Y - X*ev.B;
  arma::mat Rt = Ur.t()*R;
  arma::mat At(Rt.n_rows, Rt.n_cols);
  for(arma::uword i=0; i<Rt.n_rows; ++i){
    arma::mat Omega = d.s(i)*Vu + Ve;
    arma::vec a;
    if(!arma::solve(a, Omega, Rt.row(i).t(), arma::solve_opts::likely_sympd)){
      Rcpp::stop("Mixed model equations could not be solved");
    }
    At.row(i) = a.t();
  }
  arma::mat Py = Ur*At;
  if(d.nNull>0.0){
    Py += (R - Ur*Rt)*arma::inv_sympd(Ve);
  }
  return Py;
}

void checkInput(const arma::mat& Y, const arma::mat& X){
  if(Y.n_rows!=X.n_rows){
    Rcpp::stop("Y and X must have the same number of rows");
  }
  if(Y.n_cols<1){
    Rcpp::stop("Y must have at least one column");
  }
  if(X.n_cols<1){
    Rcpp::stop("X must have at least one column");
  }
  if(X.n_cols>=X.n_rows){
    Rcpp::stop("X must have fewer columns than rows");
  }
  if(!Y.is_finite() || !X.is_finite()){
    Rcpp::stop("Y and X must not contain missing or infinite values");
  }
}

} // namespace

//' @title Solve Multivariate RR-BLUP
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Mu+e}, where the
//' marker effects of the t traits have covariance \eqn{V_u} among traits and
//' are independent across markers, and the errors have covariance \eqn{V_e}
//' among traits and are independent across records. \eqn{V_u} and \eqn{V_e}
//' are unstructured and estimated by REML.
//'
//' @details
//' The records are rotated by the eigendecomposition of \eqn{MM'}, and the
//' traits by a second transformation that diagonalizes \eqn{V_u} and
//' \eqn{V_e} together, following \insertCite{zhou2014Efficient;textual}{AlphaSimR}.
//' The likelihood and its derivatives are then linear in the number of
//' records. When there are fewer markers than records, the eigenvectors come
//' from a thin singular value decomposition of M and the records in the null
//' space of \eqn{MM'} are summarized, so no matrix with n rows and n columns
//' is formed.
//'
//' Starting values come from univariate REML fits of each trait, by the method
//' of \insertCite{lippert2011FaST;textual}{AlphaSimR}, with the correlation
//' among least squares residuals for both covariance matrices. Five EM
//' iterations follow, and then Newton steps in the Cholesky factors of
//' \eqn{V_u} and \eqn{V_e}, whose Hessian is built from the average
//' information matrix \insertCite{gilmour1995Average}{AlphaSimR}. Working
//' with the factors keeps both matrices positive semidefinite, so an
//' estimate on the boundary of the parameter space, such as a correlation of
//' one, is reached without leaving it. An EM step is taken whenever a Newton
//' step fails to increase the likelihood. Iteration stops when the restricted
//' log-likelihood changes by less than \code{tol}.
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param maxIter maximum number of iterations
//' @param tol tolerance for convergence
//'
//' @return A list containing
//' \describe{
//'   \item{Vu}{q by q covariance matrix of the marker effects}
//'   \item{Ve}{q by q covariance matrix of the errors}
//'   \item{beta}{x by q matrix of fixed effects}
//'   \item{u}{m by q matrix of marker effects}
//'   \item{LL}{the restricted log-likelihood
//'     \eqn{-\frac{1}{2}[(N-p)\log(2\pi)+\log|V|+\log|X'V^{-1}X|+y'Py]} of
//'     \eqn{y=vec(Y)}, where N is nq, p is xq and the
//'     \eqn{\log|X'X|} term is excluded}
//'   \item{iter}{number of iterations taken}
//' }
//' A warning is given if \code{maxIter} is reached without convergence.
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X,
                         const arma::mat& M, int maxIter=1000,
                         double tol=1e-6){
  checkInput(Y, X);
  if(M.n_rows!=Y.n_rows){
    Rcpp::stop("Y and M must have the same number of rows");
  }
  arma::mat Ur;
  bool useFactor = M.n_cols<M.n_rows;
  arma::mat Kn;
  if(!useFactor){
    Kn = M*M.t();
  }
  MvData d = rotateData(Y, X, M, Kn, useFactor, Ur);
  Kn.reset();
  arma::mat Vu, Ve;
  int iter = 0;
  MvEval ev = fitMv(d, Vu, Ve, maxIter, tol, iter);
  // u = (Vu (x) M')*vec(P*y) = M'*(P*y)*Vu
  arma::mat u = M.t()*recordPy(Y, X, Ur, d, ev, Vu, Ve)*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=ev.B,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ev.LL,
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multivariate Model
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Zu+e}, where the
//' random effects have covariance \eqn{V_u \otimes K} and the errors
//' \eqn{V_e \otimes I}, with \eqn{V_u} and \eqn{V_e} the unstructured
//' covariance matrices among the q traits. They are estimated by REML.
//'
//' @details
//' The method is the one described for \code{\link{solveRRBLUPMV}}, with the
//' kernel \eqn{ZKZ'} on the records in place of \eqn{MM'}
//' \insertCite{zhou2014Efficient}{AlphaSimR}. K may be singular, and Z may
//' map several records to one level. When there are fewer levels than
//' records, K is factored as \eqn{K=FF'} and the method works with \eqn{ZF},
//' so no matrix with n rows and n columns is formed.
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iterations
//'
//' @return A list containing
//' \describe{
//'   \item{Vu}{q by q covariance matrix of the random effects}
//'   \item{Ve}{q by q covariance matrix of the errors}
//'   \item{beta}{x by q matrix of fixed effects}
//'   \item{u}{m by q matrix of random effects}
//'   \item{LL}{the restricted log-likelihood
//'     \eqn{-\frac{1}{2}[(N-p)\log(2\pi)+\log|V|+\log|X'V^{-1}X|+y'Py]} of
//'     \eqn{y=vec(Y)}, where N is nq, p is xq and the
//'     \eqn{\log|X'X|} term is excluded}
//'   \item{iter}{number of iterations taken}
//' }
//' A warning is given if \code{maxIter} is reached without convergence.
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K,
                    double tol=1e-6, int maxIter=1000){
  checkInput(Y, X);
  if(Z.n_rows!=Y.n_rows){
    Rcpp::stop("Y and Z must have the same number of rows");
  }
  if( (K.n_rows!=Z.n_cols) || (K.n_cols!=Z.n_cols) ){
    Rcpp::stop("K must be square with one row for each column of Z");
  }
  arma::mat Ks = 0.5*(K+K.t());
  arma::mat Ur, F, Kn;
  bool useFactor = Z.n_cols<Z.n_rows;
  if(useFactor){
    // K = Uk*diag(sk)*Uk', so Z*K*Z' = F*F' for F = Z*Uk*diag(sqrt(sk)).
    // Eigenvalues that are not positive carry nothing into the kernel.
    arma::uword m = Ks.n_rows;
    arma::vec sk(m);
    arma::mat Uk(m, m);
    if(eigen2(sk, Uk, Ks)!=0){
      Rcpp::stop("Decomposition of K failed");
    }
    double skMax = std::max(sk.max(), 0.0);
    arma::uvec keep = find(sk>(double(m)*arma::datum::eps*skMax));
    F = Z*Uk.cols(keep);
    F.each_row() %= sqrt(sk(keep)).t();
  }else{
    Kn = Z*Ks*Z.t();
  }
  MvData d = rotateData(Y, X, F, Kn, useFactor, Ur);
  F.reset();
  Kn.reset();
  arma::mat Vu, Ve;
  int iter = 0;
  MvEval ev = fitMv(d, Vu, Ve, maxIter, tol, iter);
  // u = (Vu (x) K*Z')*vec(P*y) = K*Z'*(P*y)*Vu
  arma::mat u = Ks*(Z.t()*recordPy(Y, X, Ur, d, ev, Vu, Ve))*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=ev.B,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ev.LL,
                            Rcpp::Named("iter")=iter);
}
