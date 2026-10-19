#include "mme.h"
#include "mmeChol.h"

// Single trait mixed models with several random terms,
//
//   y = X*beta + sum_i Z_i*u_i + e,  u_i ~ N(0, K_i*Vu_i),  e ~ N(0, I*Ve)
//
// fitted by average information REML (Gilmour, Thompson and Cullis 1995;
// Johnson and Thompson 1995). With the marginal variance
//
//   V = sum_i Vu_i*G_i + Ve*I,  G_i = Z_i*K_i*Z_i'
//
// and P = V^-1 - V^-1*X*(X'*V^-1*X)^-1*X'*V^-1, the derivative of the REML
// log-likelihood with respect to a variance component theta_i whose
// covariance matrix is G_i is
//
//   dl/dtheta_i = 0.5*(y'*P*G_i*P*y - tr(P*G_i))
//
// and the average of the observed and expected information is
//
//   AI_ij = 0.5*y'*P*G_i*P*G_j*P*y
//
// The residual is the term with G = I. Newton's method with AI in place of
// the Hessian is the update.
//
// The same quantities can be reached in two ways, and the cheaper is used.
//
// On the record side everything is done with V itself, which has the
// dimension of the records whatever the kernels are. Each G_i is held as a
// product L_i*R_i' of two thin matrices, L_i = Z_i*K_i and R_i = Z_i, or
// L_i = R_i = M_i in the RR-BLUP form. When those have fewer columns than
// there are records the n by n matrix G_i is never formed: G_i*v is
// L_i*(R_i'*v) and tr(P*G_i) = tr((P*L_i)*R_i') is the sum of the
// elementwise product of P*L_i and R_i. Otherwise forming G_i once is cheaper
// than repeating those products in every iteration, and it is kept.
//
// On the coefficient side everything is done with Henderson's mixed model
// equations. Every term is written as G_i = W_i*W_i' with effects
// a_i ~ N(0, I*Vu_i): W_i = M_i in the RR-BLUP form, and W_i = Z_i*F_i for a
// factor K_i = F_i*F_i' otherwise, with u_i = F_i*a_i. For the design
// T = [X W_1 ... W_k] the coefficient matrix is
//
//   C = T'*T/Ve + blockdiag(0, I/Vu_1, ..., I/Vu_k)
//
// and with s the solution of C*s = T'*y/Ve and e = y - T*s,
//
//   log|V| + log|X'*V^-1*X| = n*log(Ve) + sum_i r_i*log(Vu_i) + log|C|
//   P*y = e/Ve
//   tr(P*G_i) = r_i/Vu_i - tr(C^ii)/Vu_i^2
//   tr(P) = n/Ve - tr(C^-1*T'*T)/Ve^2
//   P*w = (w - T*C^-1*T'*w/Ve)/Ve
//
// where r_i is the number of columns of W_i and C^ii the diagonal block of
// C^-1 for term i. The average information comes from the k+1 working
// variates G_i*P*y and P*y by the last identity, so V and P are never formed.
// The likelihood, scores and information are the same on both sides; only
// the cost differs, and useCoefficientSide in mme.h chooses between them.

namespace {

// Fixed tolerance for solveRRBLUPMK, which takes no tol argument. It is
// applied the same way as tol in solveMKM.
const double rrblupMKTol = 1e-6;

// Most times the step length can be halved before falling back to an EM step
const int maxHalving = 10;

enum class Side { record, coefficient };

// A contiguous block of columns of T that is in the equations, and where its
// coefficients start in the solution
struct Segment {
  arma::uword start; // first column of T
  arma::uword len;
  arma::uword pos;   // first position in the solution
};

// Everything the iterations need at one value of the variance components
struct RemlPoint {
  arma::vec theta; // Vu_1..Vu_k, then Ve
  double LL = 0.0;
  bool ok = false;
  // Record side
  arma::mat Lv;    // lower Cholesky factor of V
  arma::vec Py;
  arma::mat VinvX;
  arma::mat Ainv;  // (X'*V^-1*X)^-1
  // Coefficient side. A term whose variance is zero has no effects, so it is
  // left out of the equations rather than given an infinite precision.
  std::vector<Segment> segs; // the fixed effects, then each term in the equations
  std::vector<long> termSeg; // segment of each term, or -1 if left out
  arma::uvec idx;            // columns of T in the equations
  arma::mat Lc;    // lower Cholesky factor of C
  arma::vec sol;
  arma::vec e;               // y - T*s
};

class MultiKernelREML {
public:
  // On the record side G_i = L_i*R_i'. On the coefficient side L_i = R_i = W_i
  // and the effects of the term are F_i times those of W_i, or those of W_i
  // themselves when F_i is null.
  MultiKernelREML(const arma::vec& y_, const arma::mat& X, Side side_,
                  const std::vector<const arma::mat*>& L_,
                  const std::vector<const arma::mat*>& R_,
                  const std::vector<const arma::mat*>& F_)
    : side(side_), y(y_), L(L_), R(R_), F(F_){
    k = L.size();
    setFixed(X);
    meanDiag.set_size(k);
    rankBound.set_size(k);
    confounded.zeros(k);
    for(arma::uword i=0; i<k; ++i){
      const arma::mat& Li = *L[i];
      const arma::mat& Ri = *R[i];
      meanDiag(i) = arma::accu(Li%Ri)/double(n);
      rankBound(i) = double(std::min(n, Li.n_cols));
      // A term whose G_i lies entirely in the span of X is absorbed by the
      // fixed effects: the likelihood does not depend on its variance, which
      // leaves the average information singular. The test is the trace of
      // G_i after projecting out X, against the rounding in forming it.
      arma::mat Lp = Li - Xr*(Xr.t()*Li);
      arma::mat Rp = Ri - Xr*(Xr.t()*Ri);
      double freeTrace = arma::accu(Lp%Rp);
      double roundTrace = arma::accu(arma::abs(Li)%arma::abs(Ri));
      confounded(i) = !(freeTrace > 100.0*double(n)*arma::datum::eps*roundTrace);
    }
    if(side==Side::record){
      G.set_size(k);
      dense.assign(k, false);
      for(arma::uword i=0; i<k; ++i){
        const arma::mat& Li = *L[i];
        const arma::mat& Ri = *R[i];
        if(Li.n_cols>=n){
          dense[i] = true;
          G(i) = Li*Ri.t();
          G(i) = 0.5*(G(i)+G(i).t()); // rounding only
        }
      }
    }else{
      // T = [Xr W_1 ... W_k], with T'*T and T'*y formed once
      termStart.set_size(k);
      termLen.set_size(k);
      arma::uword c = p;
      for(arma::uword i=0; i<k; ++i){
        termStart(i) = c;
        termLen(i) = L[i]->n_cols;
        c += termLen(i);
      }
      T.set_size(n, c);
      T.cols(0, p-1) = Xr;
      for(arma::uword i=0; i<k; ++i){
        if(termLen(i)>0){
          T.cols(termStart(i), termStart(i)+termLen(i)-1) = *L[i];
        }
      }
      TtT = T.t()*T;
      TtT = 0.5*(TtT+TtT.t()); // rounding only
      Tty = T.t()*y;
    }
  }

  // The record side from each G_i itself, for a caller that has already
  // formed them and needs only the variance components, so that no thin
  // factor of G_i has to be found first. The matrices are moved out of Glist.
  // rank holds an upper bound for the rank of each G_i, which the EM step
  // uses. The BLUPs need the designs, so blup() cannot be used.
  //
  // The quantities the other constructor takes from L_i and R_i follow from
  // G_i: its mean diagonal is tr(G_i)/n, and the trace left once X is
  // projected out is tr(G_i)-tr(Xr'*G_i*Xr), with Xr orthonormal.
  MultiKernelREML(const arma::vec& y_, const arma::mat& X,
                  arma::field<arma::mat>& Glist, const arma::vec& rank)
    : side(Side::record), y(y_){
    k = Glist.n_elem;
    setFixed(X);
    L.assign(k, nullptr);
    R.assign(k, nullptr);
    F.assign(k, nullptr);
    meanDiag.set_size(k);
    rankBound = rank;
    confounded.zeros(k);
    G.set_size(k);
    dense.assign(k, true);
    for(arma::uword i=0; i<k; ++i){
      if( (Glist(i).n_rows!=n) || (Glist(i).n_cols!=n) ){
        Rcpp::stop("each G must be square with one row per record");
      }
      G(i) = std::move(Glist(i));
      G(i) = 0.5*(G(i)+G(i).t()); // rounding only
      double trG = arma::trace(G(i));
      meanDiag(i) = trG/double(n);
      double freeTrace = trG - arma::accu(Xr%(G(i)*Xr));
      double roundTrace = arma::accu(arma::abs(G(i).diag()));
      confounded(i) = !(freeTrace > 100.0*double(n)*arma::datum::eps*roundTrace);
    }
  }

  // Number of records, terms and the rank of X
  arma::uword n, k, p;
  FixedDesign fd;
  Side side;

  // Starting values. The residual variance after fitting the fixed effects
  // by least squares is split in half between the residual and the random
  // terms, and the random half is shared equally among the terms on the
  // scale of the records: a term contributes Vu_i times the mean diagonal
  // of G_i to the variance of a record. This is a neutral guess that keeps
  // every component strictly positive, which an iteration needs to be able
  // to move a component in either direction. A term whose G_i has a zero
  // diagonal explains nothing and is held at zero throughout, and so is a
  // term absorbed by the fixed effects, whose variance cannot be estimated.
  arma::vec startValues(){
    arma::vec b = arma::solve(Xr, y);
    arma::vec e = y - Xr*b;
    vy = arma::dot(e,e)/double(n-p);
    if(!(vy>0.0) || !std::isfinite(vy)){
      vy = 1.0;
    }
    veFloor = 1e-8*vy;
    estimable.set_size(k);
    double nEst = 0.0;
    for(arma::uword i=0; i<k; ++i){
      estimable(i) = (meanDiag(i)>0.0) && std::isfinite(meanDiag(i)) &&
        (confounded(i)==0);
      if(estimable(i)) nEst += 1.0;
    }
    arma::vec theta(k+1, arma::fill::zeros);
    for(arma::uword i=0; i<k; ++i){
      if(estimable(i)){
        theta(i) = 0.5*vy/(nEst*meanDiag(i));
      }
    }
    theta(k) = (nEst>0.0) ? 0.5*vy : vy;
    return theta;
  }

  // Fills in the REML log-likelihood at pt.theta
  //
  //   l = -0.5*[(n-p)*log(2*pi) + log|V| + log|X'*V^-1*X| + y'*P*y]
  //
  // together with the Cholesky factor the derivatives start from. Returns
  // false if V, or the matrix that stands in for it, is not positive
  // definite. The likelihood needs only the factor, and an inverse costs
  // twice as much again, so the inverse is left to derivatives(). Most
  // evaluations are trial points in a line search, and only an accepted
  // point has its derivatives taken.
  bool evaluate(RemlPoint& pt) const{
    if(side==Side::record){
      return evaluateRecord(pt);
    }
    return evaluateCoefficient(pt);
  }

  // Score vector and average information matrix at an evaluated point
  void derivatives(const RemlPoint& pt, arma::vec& score,
                   arma::mat& AI) const{
    if(side==Side::record){
      derivativesRecord(pt, score, AI);
    }else{
      derivativesCoefficient(pt, score, AI);
    }
  }

  // Keeps every genetic component non-negative and the residual above a
  // small floor, so that V stays positive definite
  void project(arma::vec& theta) const{
    for(arma::uword i=0; i<k; ++i){
      if(!estimable(i) || !(theta(i)>0.0)){
        theta(i) = 0.0;
      }
    }
    if(!(theta(k)>veFloor)){
      theta(k) = veFloor;
    }
  }

  // An EM-REML step written in terms of the score. For a term whose kernel
  // has rank r the EM update of Henderson's equations is
  //
  //   Vu_new = Vu + Vu^2/r*(y'*P*G*P*y - tr(P*G))
  //
  // which cannot go negative, because Vu*tr(P*G) is at most r. Using an
  // upper bound for r, the smaller dimension of the term's design, keeps
  // that property and only shortens the step. EM is slow but rarely fails to
  // increase the likelihood, so it is the fallback when the Newton step
  // does not.
  arma::vec emStep(const RemlPoint& pt, const arma::vec& score) const{
    arma::vec theta = pt.theta;
    for(arma::uword i=0; i<k; ++i){
      if(estimable(i)){
        theta(i) += theta(i)*theta(i)/rankBound(i)*2.0*score(i);
      }
    }
    theta(k) += theta(k)*theta(k)/double(n)*2.0*score(k);
    project(theta);
    return theta;
  }

  // A variance for term i on the scale of the data: the one at which the term
  // alone would account for the residual variance of a record
  double typicalVu(arma::uword i) const{
    return vy/meanDiag(i);
  }

  // First-order optimality, on the boundary and off it. The gain from moving
  // a component to Vu is, to first order, score*(Vu - current). A component
  // held at zero is optimal only if moving it up to a typical variance would
  // not raise the likelihood, and one away from zero only if neither doubling
  // it, nor moving it to a typical variance, nor taking it to zero would.
  // Each gain must be below a small multiple of the tolerance. Checking the
  // free components matters as much as the boundary: a run of short steps
  // can satisfy the small-change test well away from the optimum.
  bool stationary(const RemlPoint& pt, const arma::vec& score,
                  double tol) const{
    for(arma::uword i=0; i<k; ++i){
      if(!estimable(i)) continue;
      if(pt.theta(i)==0.0){
        if(score(i)*typicalVu(i) > 10.0*tol) return false;
      }else{
        double s = std::max(pt.theta(i), typicalVu(i));
        if(std::abs(score(i))*s > 10.0*tol) return false;
      }
    }
    if(pt.theta(k)>veFloor){
      double s = std::max(pt.theta(k), vy);
      if(std::abs(score(k))*s > 10.0*tol) return false;
    }else if(score(k)*vy > 10.0*tol){
      return false;
    }
    return true;
  }

  // Projected gradient ascent, the fallback that can move a component off
  // zero. Each component moves by its score times the square of its scale,
  // which is its current value or a typical variance at zero, and the step is
  // halved until the likelihood strictly increases. A small enough step always
  // increases it when a score is nonzero, so this only fails at a stationary
  // point or when rounding prevents any gain.
  bool releaseStep(const RemlPoint& pt, const arma::vec& score,
                   RemlPoint& out) const{
    arma::vec d(k+1, arma::fill::zeros);
    for(arma::uword i=0; i<k; ++i){
      if(estimable(i) && (pt.theta(i)>0.0 || score(i)>0.0)){
        double s = (pt.theta(i)>0.0) ? pt.theta(i) : typicalVu(i);
        d(i) = score(i)*s*s;
      }
    }
    d(k) = score(k)*pt.theta(k)*pt.theta(k);
    double step = 1.0;
    for(int h=0; h<60; ++h){
      out.theta = pt.theta + step*d;
      project(out.theta);
      if(evaluate(out) && (out.LL>pt.LL)){
        return true;
      }
      step *= 0.5;
    }
    return false;
  }

  // Size of a change in the components relative to the variance of a
  // record. Each component is weighted by the mean diagonal of its G, so a
  // change is judged by how much variance it moves rather than by the
  // arbitrary scale of the marker codes, and a component near zero does not
  // make the criterion meaningless.
  double relativeChange(const arma::vec& a, const arma::vec& b) const{
    arma::vec w(k+1);
    w.head(k) = meanDiag;
    w(k) = 1.0;
    for(arma::uword i=0; i<k; ++i){
      if(!estimable(i)) w(i) = 0.0;
    }
    return arma::accu(arma::abs(a-b)%w)/arma::accu(arma::abs(b)%w);
  }

  // Runs the iterations from the starting values and leaves the final point
  // in pt. Returns the number of iterations used.
  int fit(RemlPoint& pt, int maxIter, double tol, const char* name){
    pt.theta = startValues();
    if(!evaluate(pt)){
      Rcpp::stop("%s: the variance matrix is not positive definite at the "
                   "starting values", name);
    }
    int iter = 0;
    bool converged = false;
    arma::vec score;
    arma::mat AI;
    while(iter<maxIter){
      ++iter;
      derivatives(pt, score, AI);

      // A component on the boundary stays there when the likelihood would
      // fall by moving it up, and is dropped from the Newton system. This is
      // what lets a dominance or SCA variance settle at exactly zero instead
      // of bouncing off it.
      std::vector<arma::uword> freeIdx;
      for(arma::uword i=0; i<k; ++i){
        if(estimable(i) && (pt.theta(i)>0.0 || score(i)>0.0)){
          freeIdx.push_back(i);
        }
      }
      freeIdx.push_back(k);
      arma::uvec f(freeIdx);
      arma::vec delta(k+1, arma::fill::zeros);
      arma::mat Lai;
      bool newtonOk = arma::chol(Lai, AI.submat(f,f), "lower");
      if(newtonOk){
        arma::vec sf = score.elem(f);
        arma::vec df = arma::solve(arma::trimatu(Lai.t()),
                                   arma::solve(arma::trimatl(Lai), sf));
        delta.elem(f) = df;
        newtonOk = delta.is_finite();
      }

      // Accept a step that does not lower the likelihood, allowing for
      // rounding in a likelihood that is flat at the optimum
      double slack = 1e-10*(1.0+std::abs(pt.LL));
      RemlPoint trial;
      bool accepted = false;
      if(newtonOk){
        double step = 1.0;
        for(int h=0; h<=maxHalving; ++h){
          trial.theta = pt.theta + step*delta;
          project(trial.theta);
          if(evaluate(trial) && trial.LL>=pt.LL-slack){
            accepted = true;
            break;
          }
          step *= 0.5;
        }
      }
      if(!accepted){
        trial.theta = emStep(pt, score);
        accepted = evaluate(trial) && trial.LL>=pt.LL-slack;
      }
      if(accepted){
        double dLL = trial.LL - pt.LL;
        double dTheta = relativeChange(trial.theta, pt.theta);
        pt = trial;
        if(!(std::abs(dLL)<tol && dTheta<tol)){
          continue;
        }
        // A small change is only convergence at a stationary point, where no
        // component, at zero or not, could raise the likelihood by moving
        derivatives(pt, score, AI);
        if(stationary(pt, score, tol)){
          converged = true;
          break;
        }
      }else if(stationary(pt, score, tol)){
        // No step improves the likelihood and the scores say none would, so
        // this is a stationary point to within the precision of the
        // arithmetic
        converged = true;
        break;
      }
      // The scores show the likelihood can still rise but the last step made
      // no real progress. Typically a component at zero has a positive score,
      // which the EM step cannot act on because it scales the change by the
      // square of the component. A projected gradient step can.
      RemlPoint released;
      if(!releaseStep(pt, score, released)){
        // Not even a step of 2^-60 along the scores raises the likelihood,
        // so what the scores show is rounding: the point is stationary to
        // the precision of the arithmetic
        converged = true;
        break;
      }
      pt = released;
    }
    if(!converged){
      Rf_warning("%s: did not converge in %d iterations", name, maxIter);
    }
    return iter;
  }

  // Generalized least squares estimates, with zeros for any columns of X
  // that were dropped for being linearly dependent on the others
  arma::mat beta(const RemlPoint& pt, arma::uword nColX) const{
    arma::vec br;
    if(side==Side::record){
      br = pt.Ainv*(pt.VinvX.t()*y);
    }else{
      br = pt.sol.head(p);
    }
    return fixedEffectsFromBasis(fd, br+b0, nColX);
  }

  // BLUP of the term's effects. On the record side this is
  // u_i = Vu_i*K_i*Z_i'*P*y = Vu_i*L_i'*P*y. On the coefficient side it is
  // read from the solution, mapped through F_i where the term has one.
  arma::mat blup(const RemlPoint& pt, arma::uword i) const{
    if(side==Side::record){
      if(L[i]==nullptr){
        Rcpp::stop("the BLUPs need the designs, which were not given");
      }
      arma::mat u = pt.theta(i)*(L[i]->t()*pt.Py);
      return u;
    }
    arma::vec a(termLen(i), arma::fill::zeros);
    if(pt.termSeg[i]>=0){
      const Segment& s = pt.segs[pt.termSeg[i]];
      a = pt.sol.subvec(s.pos, s.pos+s.len-1);
    }
    if(F[i]!=nullptr){
      arma::mat u = (*F[i])*a;
      return u;
    }
    arma::mat u = a;
    return u;
  }

private:
  arma::vec y;  // the response with its least squares fit taken off
  arma::vec b0; // that least squares fit
  arma::mat Xr;
  std::vector<const arma::mat*> L, R, F;
  arma::vec meanDiag, rankBound;
  arma::uvec estimable, confounded;
  double vy = 1.0, veFloor = 0.0;
  // Record side
  arma::field<arma::mat> G;
  std::vector<bool> dense;
  // Coefficient side
  arma::mat T, TtT;
  arma::vec Tty;
  arma::uvec termStart, termLen;

  // Sets up the fixed effects and takes their least squares fit off y.
  // Columns of X that depend on the others are dropped, and the fit uses Xr,
  // an orthonormal basis for the span of the rest, so a badly scaled X costs
  // no precision; see fixedDesign. The least squares fit is added back to the
  // fixed effects by beta(), and taking it off keeps a large mean from
  // swamping the sums of squares.
  void setFixed(const arma::mat& X){
    n = y.n_elem;
    fd = fixedDesign(X);
    if(fd.keep.n_elem==0){
      Rcpp::stop("X must have at least one nonzero column");
    }
    Xr = fd.Q;
    p = Xr.n_cols;
    if(p>=n){
      Rcpp::stop("the rank of X must be smaller than the number of records");
    }
    arma::mat B0;
    y = removeLeastSquares(Xr, y, B0).col(0);
    b0 = B0.col(0);
  }

  bool evaluateRecord(RemlPoint& pt) const{
    pt.ok = false;
    arma::mat V(n, n, arma::fill::zeros);
    V.diag().fill(pt.theta(k));
    for(arma::uword i=0; i<k; ++i){
      if(pt.theta(i)==0.0) continue;
      if(dense[i]){
        V += pt.theta(i)*G(i);
      }else{
        V += pt.theta(i)*((*L[i])*(R[i]->t()));
      }
    }
    V = 0.5*(V+V.t()); // rounding only
    if(!arma::chol(pt.Lv, V, "lower")){
      return false;
    }
    V.reset();
    double logDetV = 2.0*arma::accu(arma::log(pt.Lv.diag()));
    // With V = Lv*Lv', Q = Lv^-1*X and z = Lv^-1*y give X'*V^-1*X = Q'*Q,
    // and r = z - Q*(Q'*Q)^-1*Q'*z gives P*y = Lv'^-1*r and y'*P*y = r'*r
    const arma::solve_opts::opts fast = arma::solve_opts::fast;
    arma::mat Q = arma::solve(arma::trimatl(pt.Lv), Xr, fast);
    arma::vec z = arma::solve(arma::trimatl(pt.Lv), y, fast);
    arma::mat A = Q.t()*Q;
    A = 0.5*(A+A.t()); // rounding only
    double logDetA;
    if(!invAndLogDet(pt.Ainv, logDetA, A)){
      return false;
    }
    arma::vec r = z - Q*(pt.Ainv*(Q.t()*z));
    pt.Py = arma::solve(arma::trimatu(pt.Lv.t()), r, fast);
    pt.VinvX = arma::solve(arma::trimatu(pt.Lv.t()), Q, fast);
    double yPy = arma::dot(r, r);
    pt.LL = -0.5*(double(n-p)*std::log(2.0*arma::datum::pi) + logDetV +
      logDetA + yPy);
    pt.ok = std::isfinite(pt.LL);
    return pt.ok;
  }

  void derivativesRecord(const RemlPoint& pt, arma::vec& score,
                         arma::mat& AI) const{
    score.set_size(k+1);
    // P = V^-1 - V^-1*X*(X'*V^-1*X)^-1*X'*V^-1, with V^-1 = Lv'^-1*Lv^-1
    arma::mat P;
    {
      arma::mat Linv;
      if(!arma::inv(Linv, arma::trimatl(pt.Lv))){
        Rcpp::stop("could not invert the Cholesky factor of V");
      }
      P = Linv.t()*Linv;
    }
    P -= pt.VinvX*pt.Ainv*pt.VinvX.t();
    P = 0.5*(P+P.t()); // rounding only
    arma::mat W(n, k+1);
    for(arma::uword i=0; i<k; ++i){
      double trPG;
      if(dense[i]){
        W.col(i) = G(i)*pt.Py;
        trPG = arma::accu(P%G(i));
      }else{
        W.col(i) = (*L[i])*(R[i]->t()*pt.Py);
        trPG = arma::accu((P*(*L[i]))%(*R[i]));
      }
      score(i) = 0.5*(arma::dot(pt.Py, W.col(i)) - trPG);
    }
    W.col(k) = pt.Py;
    score(k) = 0.5*(arma::dot(pt.Py, pt.Py) - arma::trace(P));
    AI = 0.5*(W.t()*(P*W));
    AI = 0.5*(AI+AI.t()); // rounding only
  }

  // T_a'*M, where T_a holds the columns of T in the equations
  arma::mat crossT(const RemlPoint& pt, const arma::mat& M) const{
    arma::mat out(pt.idx.n_elem, M.n_cols);
    for(const Segment& s : pt.segs){
      out.rows(s.pos, s.pos+s.len-1) =
        T.cols(s.start, s.start+s.len-1).t()*M;
    }
    return out;
  }

  // T_a*S
  arma::mat multT(const RemlPoint& pt, const arma::mat& S) const{
    arma::mat out(n, S.n_cols, arma::fill::zeros);
    for(const Segment& s : pt.segs){
      out += T.cols(s.start, s.start+s.len-1)*S.rows(s.pos, s.pos+s.len-1);
    }
    return out;
  }

  bool evaluateCoefficient(RemlPoint& pt) const{
    pt.ok = false;
    const double ve = pt.theta(k);
    pt.segs.clear();
    pt.termSeg.assign(k, -1);
    pt.segs.push_back(Segment{0, p, 0});
    arma::uword pos = p;
    double logDetG = 0.0;
    for(arma::uword i=0; i<k; ++i){
      if( (pt.theta(i)>0.0) && (termLen(i)>0) ){
        pt.termSeg[i] = long(pt.segs.size());
        pt.segs.push_back(Segment{termStart(i), termLen(i), pos});
        pos += termLen(i);
        logDetG += double(termLen(i))*std::log(pt.theta(i));
      }
    }
    pt.idx.set_size(pos);
    for(const Segment& s : pt.segs){
      for(arma::uword j=0; j<s.len; ++j){
        pt.idx(s.pos+j) = s.start+j;
      }
    }
    arma::mat C = TtT.submat(pt.idx, pt.idx)/ve;
    for(arma::uword i=0; i<k; ++i){
      if(pt.termSeg[i]>=0){
        const Segment& s = pt.segs[pt.termSeg[i]];
        for(arma::uword j=s.pos; j<s.pos+s.len; ++j){
          C(j,j) += 1.0/pt.theta(i);
        }
      }
    }
    if(!arma::chol(pt.Lc, C, "lower")){
      return false;
    }
    C.reset();
    double logDetC = 2.0*arma::accu(arma::log(pt.Lc.diag()));
    arma::vec rhs = Tty.elem(pt.idx)/ve;
    const arma::solve_opts::opts fast = arma::solve_opts::fast;
    pt.sol = arma::solve(arma::trimatu(pt.Lc.t()),
                         arma::solve(arma::trimatl(pt.Lc), rhs, fast), fast);
    pt.e = y - multT(pt, pt.sol);
    double yPy = arma::dot(y, pt.e)/ve;
    pt.LL = -0.5*(double(n-p)*std::log(2.0*arma::datum::pi) +
      double(n)*std::log(ve) + logDetG + logDetC + yPy);
    pt.ok = std::isfinite(pt.LL);
    return pt.ok;
  }

  void derivativesCoefficient(const RemlPoint& pt, arma::vec& score,
                              arma::mat& AI) const{
    const double ve = pt.theta(k);
    score.set_size(k+1);
    // C^-1 = Lc'^-1*Lc^-1 is never formed. The traces need only its diagonal
    // blocks, and the trace of a block is the sum of squares of the matching
    // columns of Lc^-1, which saves the third of the cost of an inverse that
    // goes on multiplying Lc^-1 by its transpose.
    arma::mat Linv;
    if(!arma::inv(Linv, arma::trimatl(pt.Lc))){
      Rcpp::stop("could not invert the Cholesky factor of C");
    }
    // Sum over the terms in the equations of tr(C^ii)/Vu_i
    double trCD = 0.0;
    // Working variates G_i*P*y, then P*y
    arma::mat W(n, k+1, arma::fill::zeros);
    for(arma::uword i=0; i<k; ++i){
      if(termLen(i)==0){
        score(i) = 0.0;
        continue;
      }
      arma::uword a = termStart(i), b = termStart(i)+termLen(i)-1;
      arma::vec Wte = T.cols(a, b).t()*pt.e;
      W.col(i) = T.cols(a, b)*Wte/ve;
      double quad = arma::dot(Wte, Wte)/(ve*ve); // y'*P*G_i*P*y
      double trPG;
      if(pt.termSeg[i]>=0){
        const Segment& s = pt.segs[pt.termSeg[i]];
        double th = pt.theta(i);
        double trCii = arma::accu(arma::square(Linv.cols(s.pos,
                                                         s.pos+s.len-1)));
        trCD += trCii/th;
        trPG = double(s.len)/th - trCii/(th*th);
      }else{
        // A term left out of the equations: tr(P*G_i) is
        // tr(W_i'*W_i)/Ve - tr(B*C^-1*B')/Ve^2 for B = W_i'*T_a
        arma::uvec rows = arma::regspace<arma::uvec>(a, b);
        arma::mat B = TtT.submat(rows, pt.idx);
        trPG = arma::trace(TtT.submat(a, a, b, b))/ve -
          arma::accu(arma::square(Linv*B.t()))/(ve*ve);
      }
      score(i) = 0.5*(quad - trPG);
    }
    W.col(k) = pt.e/ve;
    // T_a'*T_a = Ve*(C - D), where D holds 1/Vu_i on the diagonal of each
    // term in the equations, so tr(T_a'*T_a*C^-1) = Ve*(c - tr(D*C^-1)) and
    // tr(P) = (n - c + tr(D*C^-1))/Ve
    double trP = (double(n) - double(pt.idx.n_elem) + trCD)/ve;
    score(k) = 0.5*(arma::dot(pt.e, pt.e)/(ve*ve) - trP);
    // P*W = (W - T_a*C^-1*T_a'*W/Ve)/Ve
    arma::mat PW = (W - multT(pt, Linv.t()*(Linv*crossT(pt, W)))/ve)/ve;
    AI = 0.5*(W.t()*PW);
    AI = 0.5*(AI+AI.t()); // rounding only
  }
};

void checkResponse(const arma::mat& y, const arma::mat& X, const char* name){
  if(y.n_cols!=1){
    Rcpp::stop("%s: y must have exactly one column", name);
  }
  if(X.n_rows!=y.n_rows){
    Rcpp::stop("%s: X and y must have the same number of rows", name);
  }
  if(!y.is_finite() || !X.is_finite()){
    Rcpp::stop("%s: y and X must not contain missing or infinite values",
               name);
  }
}

void checkMlist(const arma::mat& y, const arma::field<arma::mat>& Mlist,
                const char* name){
  if(Mlist.n_elem==0){
    Rcpp::stop("%s: Mlist must contain at least one matrix", name);
  }
  for(arma::uword i=0; i<Mlist.n_elem; ++i){
    if(Mlist(i).n_rows!=y.n_rows){
      Rcpp::stop("%s: every matrix in Mlist must have one row per record",
                 name);
    }
    if(!Mlist(i).is_finite()){
      Rcpp::stop("%s: Mlist must not contain missing or infinite values",
                 name);
    }
  }
}

void checkZKlist(const arma::mat& y, const arma::field<arma::mat>& Zlist,
                 const arma::field<arma::mat>& Klist, const char* name){
  if(Zlist.n_elem==0){
    Rcpp::stop("%s: Zlist must contain at least one matrix", name);
  }
  if(Klist.n_elem!=Zlist.n_elem){
    Rcpp::stop("%s: Zlist and Klist must have the same length", name);
  }
  for(arma::uword i=0; i<Zlist.n_elem; ++i){
    const arma::mat& Z = Zlist(i);
    const arma::mat& K = Klist(i);
    if(Z.n_rows!=y.n_rows){
      Rcpp::stop("%s: every matrix in Zlist must have one row per record",
                 name);
    }
    if(K.n_rows!=Z.n_cols || K.n_cols!=Z.n_cols){
      Rcpp::stop("%s: each K must be square with one row per column of Z",
                 name);
    }
    if(!Z.is_finite() || !K.is_finite()){
      Rcpp::stop("%s: Zlist and Klist must not contain missing or infinite "
                   "values", name);
    }
  }
}

// Fits either form of the model. In the RR-BLUP form (rrForm) the matrices
// in Zlist are the marker matrices and Klist is not used. side is 0 to
// choose the cheaper side, 1 for the record side and 2 for the coefficient
// side. The inputs must already have been checked.
Rcpp::List fitMultiKernel(const arma::mat& y, const arma::mat& X,
                          const arma::field<arma::mat>& Zlist,
                          const arma::field<arma::mat>& Klist, bool rrForm,
                          int maxIter, double tol, int sideRequest,
                          const char* name, bool reportSide){
  arma::uword n = y.n_rows;
  arma::uword k = Zlist.n_elem;

  // The number of coefficients is not known until the covariance matrices
  // are factored, so the side is chosen on the number of columns, which is
  // an upper bound for it. Factoring a covariance matrix adds the cube of
  // its size to the coefficient side.
  Side side;
  if(sideRequest==1){
    side = Side::record;
  }else if(sideRequest==2){
    side = Side::coefficient;
  }else{
    arma::uword c = X.n_cols;
    double extra = 0.0;
    for(arma::uword i=0; i<k; ++i){
      double m = double(Zlist(i).n_cols);
      c += Zlist(i).n_cols;
      if(!rrForm) extra += m*m*m;
    }
    side = useCoefficientSide(n, c, k, extra) ? Side::coefficient :
      Side::record;
  }

  arma::field<arma::mat> D(k), Fs(k);
  std::vector<const arma::mat*> L(k), R(k), F(k, nullptr);
  for(arma::uword i=0; i<k; ++i){
    if(rrForm){
      L[i] = &Zlist(i);
      R[i] = &Zlist(i);
      continue;
    }
    const arma::mat& Z = Zlist(i);
    // K_i is symmetrized first so that Z_i*K_i*Z_i' is symmetric
    arma::mat K = 0.5*(Klist(i)+Klist(i).t());
    if(side==Side::record){
      // K_i is checked on its own, as on the coefficient side. A check of V
      // alone would let an indefinite K_i through whenever the residual
      // variance happened to keep V positive definite. Negative eigenvalues
      // at rounding level are set to zero, as the coefficient side does by
      // factoring over the positive ones, so both sides fit the same model.
      K = psdKernel(K, std::string(name) +
                    ": each K must be positive semidefinite");
      // Z_i*K_i is formed once, and the BLUP K_i*Z_i'*P*y can be read off
      // as (Z_i*K_i)'*P*y
      D(i) = Z*K;
      L[i] = &D(i);
      R[i] = &Zlist(i);
    }else{
      // K_i = F_i*F_i' with F_i = Q*diag(sqrt(lambda)) over the positive
      // eigenvalues, so a singular K_i is handled by having fewer columns.
      // Every positive eigenvalue is kept, however small next to the
      // largest: Z_i can scale a direction back up, so whether it matters is
      // only known after Z_i is applied, and dropping it here would fit a
      // different model from the record side.
      arma::uword m = K.n_rows;
      arma::vec lambda(m);
      arma::mat Q(m, m);
      if(m>0){
        if(eigen2(lambda, Q, K)!=0){
          Rcpp::stop("%s: eigendecomposition of K failed", name);
        }
        double maxLambda = arma::abs(lambda).max();
        if(lambda.min() < -std::sqrt(arma::datum::eps)*maxLambda){
          Rcpp::stop("%s: each K must be positive semidefinite", name);
        }
        arma::uvec keep = find(lambda>0.0);
        Fs(i) = Q.cols(keep);
        Fs(i).each_row() %= sqrt(lambda(keep)).t();
      }else{
        Fs(i).set_size(0, 0);
      }
      D(i) = Z*Fs(i);
      L[i] = &D(i);
      R[i] = &D(i);
      F[i] = &Fs(i);
    }
  }

  arma::vec yv = y.col(0);
  MultiKernelREML reml(yv, X, side, L, R, F);
  RemlPoint pt;
  int iter = reml.fit(pt, maxIter, tol, name);

  arma::field<arma::mat> u(k);
  for(arma::uword i=0; i<k; ++i){
    u(i) = reml.blup(pt, i);
  }
  arma::mat Vu = pt.theta.head(k);
  arma::mat Ve(1, 1);
  Ve(0,0) = pt.theta(k);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                                      Rcpp::Named("Ve")=Ve,
                                      Rcpp::Named("beta")=reml.beta(pt, X.n_cols),
                                      Rcpp::Named("u")=u,
                                      Rcpp::Named("LL")=pt.LL-reml.fd.logDetR,
                                      Rcpp::Named("iter")=iter);
  if(reportSide){
    out.push_back((side==Side::record) ? 1 : 2, "side");
  }
  return out;
}

// Multiple trait RR-BLUP models with several random terms,
//
//   Y = X*B + sum_i M_i*U_i + E,
//   vec(U_i) ~ N(0, Vu_i (x) I),  vec(E) ~ N(0, Ve (x) I),
//
// for t traits recorded on every one of n records, with Vu_i and Ve
// unstructured t by t covariance matrices. The variance components are
// estimated by REML, by the same scheme as the single kernel multiple trait
// solver in mmeMultiTrait.cpp: a few EM steps, then Newton steps in the lower
// Cholesky factors of the covariance matrices, whose Hessian is the average
// information matrix of Gilmour et al. (1995) carried through the Jacobian
// of V = L*L' together with the exact second order term of that map.
// Working with the factors keeps every matrix positive semidefinite, so an
// estimate on the boundary of the parameter space, such as a genetic
// correlation of one or a term with no variance, is reached without leaving
// it.
//
// The rotation that separates the traits of a single kernel model has no
// counterpart here, because no transformation diagonalizes two kernels and
// the identity at once. The likelihood is therefore reached through one of
// the two routes of the single trait solver above, whichever is cheaper.
//
// On the record side the variance of vec(Y),
//
//   V = sum_i Vu_i (x) G_i + Ve (x) I,  G_i = M_i*M_i',
//
// is formed and factored, with nt rows.
//
// On the coefficient side the effects of each term are whitened, U_i =
// A_i*L_i' for Vu_i = L_i*L_i', so that vec(A_i) ~ N(0, I). The design of
// all the effects is W = [I (x) X, L_1 (x) M_1, ..., L_k (x) M_k], and with
// R = Ve (x) I Henderson's equations have the coefficient matrix
//
//   C = W'*R^-1*W + blockdiag(0, I, ..., I),
//
// whose blocks are (L_r'*Ve^-1*L_s) (x) T_r'*T_s for the parts T_r of
// T = [X M_1 ... M_k], with L_r = I for the fixed effects. T'*T is formed
// once. With s the solution of C*s = W'*R^-1*y and E = Y - W*s,
//
//   log|V| + log|X'*V^-1*X| = n*log|Ve| + log|C|
//   P*y = vec(E*Ve^-1)
//
// The identity prior of the whitened effects stays valid when a factor is
// singular, which is why the whitened form is used: a direction with no
// genetic variance simply has no effects. The derivatives of the likelihood
// with respect to an element of Vu_i or Ve need the traces
// tr(P*(E_ab (x) G)) for the unit matrices E_ab, and with
// R^-1*(E_ab (x) G)*R^-1 = (Ve^-1*E_ab*Ve^-1) (x) G those are
//
//   tr(P*(E_ab (x) G)) = tr(Ve^-1*E_ab)*tr(G) - tr(C^-1*W'*(F (x) G)*W)
//
// for F = Ve^-1*E_ab*Ve^-1. The second trace has blocks
// (L_r'*F*L_s) (x) T_r'*G*T_s, so it is a sum over pairs of terms of small
// t by t products against traces of C^-1 that do not depend on a or b, and
// those are worked out once per kernel.
//
// The data are held as matrices with one column per trait, so the layout
// assumes every trait is recorded on every record. Supporting missing
// records would mean selecting the observed entries of vec(Y) on the record
// side and summing the residual terms of the coefficient side over the
// patterns of missing traits.

// Everything the iterations need at one value of the variance components
struct MvMkPoint {
  std::vector<arma::mat> Lf; // lower Cholesky factors, the terms then Ve
  std::vector<arma::mat> V;  // Lf*Lf'
  double LL = 0.0;
  bool ok = false;
  arma::mat A;     // P*y as an n by t matrix
  // Record side
  arma::mat Lv;    // lower Cholesky factor of V
  arma::mat VinvX;
  arma::mat Ainv;  // (X'*V^-1*X)^-1
  // Coefficient side
  arma::mat Lc;    // lower Cholesky factor of C
  arma::mat VeInv;
  arma::vec sol;
};

class MvMultiKernelREML {
public:
  MvMultiKernelREML(const arma::mat& Yin, const arma::mat& X, Side side_,
                    const arma::field<arma::mat>& Mlist)
    : side(side_){
    k = Mlist.n_elem;
    setResponses(Yin, X);
    // Each term is fitted on the scale of the mean diagonal of its G_i, so the
    // Newton steps and the eigenvalue floor on their Hessian behave the same
    // whatever the units of the markers. Vu_i is returned to those units at
    // the end, and tr(G_i) is then n for every term.
    md.zeros(k);
    rankBound.set_size(k);
    held.zeros(k);
    Ms.set_size(k);
    for(arma::uword i=0; i<k; ++i){
      const arma::mat& M = Mlist(i);
      md(i) = arma::accu(arma::square(M))/double(n);
      rankBound(i) = double(std::min(n, M.n_cols));
      if( !(md(i)>0.0) || !std::isfinite(md(i)) ){
        held(i) = 1;
        Ms(i) = M;
        continue;
      }
      Ms(i) = M/std::sqrt(md(i));
      // A term that lies entirely in the span of X is absorbed by the fixed
      // effects and its variance cannot be estimated, so it is held at zero,
      // as in the single trait solver
      arma::mat Mp = Ms(i) - Xr*(Xr.t()*Ms(i));
      double freeTrace = arma::accu(arma::square(Mp));
      double roundTrace = arma::accu(arma::square(Ms(i)));
      if(!(freeTrace > 100.0*double(n)*arma::datum::eps*roundTrace)){
        held(i) = 1;
      }
    }
    paramIndex(nT, rowIdx, colIdx);
    nPar = rowIdx.n_elem;
    if(side==Side::record){
      G.set_size(k);
      for(arma::uword i=0; i<k; ++i){
        G(i) = Ms(i)*Ms(i).t();
        G(i) = 0.5*(G(i)+G(i).t()); // rounding only
      }
      Xf = arma::kron(arma::eye(nT, nT), Xr);
      yv = arma::vectorise(Y);
    }else{
      // T = [Xr M_1 ... M_k], with T'*T, T'*Y and each T'*G_i*T formed once
      nSeg = k+1;
      segStart.set_size(nSeg);
      segLen.set_size(nSeg);
      segStart(0) = 0;
      segLen(0) = p;
      arma::uword c = p;
      for(arma::uword i=0; i<k; ++i){
        segStart(i+1) = c;
        segLen(i+1) = Ms(i).n_cols;
        c += Ms(i).n_cols;
      }
      nCoef = c;
      T.set_size(n, c);
      T.cols(0, p-1) = Xr;
      for(arma::uword i=0; i<k; ++i){
        if(segLen(i+1)>0){
          T.cols(segStart(i+1), segStart(i+1)+segLen(i+1)-1) = Ms(i);
        }
      }
      TtT = T.t()*T;
      TtT = 0.5*(TtT+TtT.t()); // rounding only
      TtY = T.t()*Y;
      GG.set_size(k);
      for(arma::uword i=0; i<k; ++i){
        if(segLen(i+1)==0){
          GG(i).zeros(c, c);
          continue;
        }
        arma::mat TtM = TtT.cols(segStart(i+1), segStart(i+1)+segLen(i+1)-1);
        GG(i) = TtM*TtM.t();
        GG(i) = 0.5*(GG(i)+GG(i).t()); // rounding only
      }
    }
  }

  // The record side from each G_i = M_i*M_i' itself, for a caller that has
  // already formed them and needs only the variance components, as for the
  // single trait solver. The matrices are moved out of Glist, and rank holds
  // an upper bound for the rank of each. The BLUPs need the marker codes, so
  // blupScaled() cannot be used. The scaling and the test for a term absorbed
  // by the fixed effects follow from traces, as tr(M_i'*M_i) = tr(G_i) and
  // tr(M_i'*(I-Xr*Xr')*M_i) = tr(G_i)-tr(Xr'*G_i*Xr).
  MvMultiKernelREML(const arma::mat& Yin, const arma::mat& X,
                    arma::field<arma::mat>& Glist, const arma::vec& rank)
    : side(Side::record){
    k = Glist.n_elem;
    setResponses(Yin, X);
    haveCodes = false;
    md.zeros(k);
    rankBound = rank;
    held.zeros(k);
    Ms.set_size(k);
    G.set_size(k);
    for(arma::uword i=0; i<k; ++i){
      if( (Glist(i).n_rows!=n) || (Glist(i).n_cols!=n) ){
        Rcpp::stop("each G must be square with one row per record");
      }
      G(i) = std::move(Glist(i));
      G(i) = 0.5*(G(i)+G(i).t()); // rounding only
      md(i) = arma::trace(G(i))/double(n);
      if( !(md(i)>0.0) || !std::isfinite(md(i)) ){
        held(i) = 1;
        continue;
      }
      G(i) /= md(i);
      double roundTrace = arma::trace(G(i));
      double freeTrace = roundTrace - arma::accu(Xr%(G(i)*Xr));
      if(!(freeTrace > 100.0*double(n)*arma::datum::eps*roundTrace)){
        held(i) = 1;
      }
    }
    paramIndex(nT, rowIdx, colIdx);
    nPar = rowIdx.n_elem;
    Xf = arma::kron(arma::eye(nT, nT), Xr);
    yv = arma::vectorise(Y);
  }

  arma::uword n, nT, k, N, p;
  FixedDesign fd;
  TraitScale ts;
  arma::mat B0;
  double veFloor = 0.0;

  // Raises any eigenvalue of Ve below the floor to it
  arma::mat floorVe(const arma::mat& Ve) const{
    arma::vec val;
    arma::mat vec;
    if(!arma::eig_sym(val, vec, 0.5*(Ve+Ve.t())) || (val.min()>=veFloor)){
      return Ve;
    }
    val.elem(find(val<veFloor)).fill(veFloor);
    arma::mat out = vec*arma::diagmat(val)*vec.t();
    return 0.5*(out+out.t());
  }

  // Whether Ve is at or above the floor, allowing for the rounding of a
  // matrix that was floored exactly
  bool aboveFloor(const arma::mat& Ve) const{
    arma::vec val;
    if(!arma::eig_sym(val, 0.5*(Ve+Ve.t()))){
      return false;
    }
    return val.min() >= 0.5*veFloor;
  }

  // Starting values. Half of each trait's residual variance after the least
  // squares fit goes to the error and half is shared equally among the terms
  // that can be estimated, as in the single trait solver, and every matrix
  // takes the correlation among the least squares residuals, which is the
  // phenotypic correlation net of the fixed effects. If that leaves a matrix
  // singular the matrices start diagonal instead.
  void startValues(MvMkPoint& pt) const{
    arma::rowvec vy = arma::sum(arma::square(Y), 0)/double(n-p);
    for(arma::uword j=0; j<nT; ++j){
      if(!(vy(j)>0.0) || !std::isfinite(vy(j))){
        vy(j) = 1.0;
      }
    }
    double nEst = double(k) - double(arma::accu(held));
    arma::mat Rr = Y.t()*Y;
    arma::vec sd = sqrt(arma::clamp(Rr.diag(), arma::datum::eps,
                                    arma::datum::inf));
    arma::mat rho = Rr/(sd*sd.t());
    rho.diag().ones();
    arma::vec ve = vy.t()*((nEst>0.0) ? 0.5 : 1.0);
    arma::vec vu = vy.t()*((nEst>0.0) ? 0.5/nEst : 0.0);
    arma::mat Ve = arma::diagmat(sqrt(ve))*rho*arma::diagmat(sqrt(ve));
    arma::mat Vu = arma::diagmat(sqrt(vu))*rho*arma::diagmat(sqrt(vu));
    arma::mat L;
    if(!arma::chol(L, Ve, "lower") || ((nEst>0.0) && !arma::chol(L, Vu, "lower"))){
      Ve = arma::diagmat(ve);
      Vu = arma::diagmat(vu);
    }
    pt.V.assign(k+1, arma::mat());
    pt.Lf.assign(k+1, arma::mat());
    for(arma::uword i=0; i<k; ++i){
      pt.V[i] = held(i) ? arma::mat(nT, nT, arma::fill::zeros) : Vu;
      pt.Lf[i] = held(i) ? arma::mat(nT, nT, arma::fill::zeros) :
        lowerChol(pt.V[i]);
    }
    pt.V[k] = floorVe(Ve);
    pt.Lf[k] = lowerChol(pt.V[k]);
    pt.V[k] = pt.Lf[k]*pt.Lf[k].t();
  }

  // Fills in the REML log-likelihood at pt.V, with the factors the derivatives
  // start from. Returns false if V, or C, is not positive definite.
  bool evaluate(MvMkPoint& pt) const{
    if(!aboveFloor(pt.V[k])){
      pt.ok = false;
      return false;
    }
    if(side==Side::record){
      return evaluateRecord(pt);
    }
    return evaluateCoefficient(pt);
  }

  // Scores and average information matrix with respect to the elements on and
  // below the diagonal of each Vu_i and of Ve, listed term by term with Ve
  // last and in the order of paramIndex within each matrix
  void derivatives(const MvMkPoint& pt, arma::vec& score, arma::mat& AI) const{
    if(side==Side::record){
      derivativesRecord(pt, score, AI);
    }else{
      derivativesCoefficient(pt, score, AI);
    }
  }

  // Runs the iterations and leaves the final point in pt. Returns the number
  // of iterations used.
  int fit(MvMkPoint& pt, int maxIter, double tol, const char* name){
    const int nEM = 5;
    const arma::uword nMat = k+1;
    const arma::uword nTot = nMat*nPar;
    // The parameters that move: those of every term that can be estimated,
    // and of Ve
    std::vector<arma::uword> freeMat, freeIdx;
    for(arma::uword g=0; g<nMat; ++g){
      if( (g<k) && held(g) ) continue;
      freeMat.push_back(g);
      for(arma::uword q=0; q<nPar; ++q){
        freeIdx.push_back(g*nPar+q);
      }
    }
    arma::uvec f(freeIdx);

    startValues(pt);
    if(!evaluate(pt)){
      Rcpp::stop("%s: the variance matrix is not positive definite at the "
                   "starting values", name);
    }
    int iter = 0;
    bool converged = false;
    arma::vec gV;
    arma::mat AI;
    for(iter=1; iter<=maxIter; ++iter){
      derivatives(pt, gV, AI);
      MvMkPoint next;
      bool accepted = false;
      if(iter>nEM){
        arma::mat J(nTot, nTot, arma::fill::zeros);
        arma::mat Kc(nTot, nTot, arma::fill::zeros);
        for(arma::uword g : freeMat){
          arma::uword a = g*nPar, b = g*nPar+nPar-1;
          J.submat(a, a, b, b) = cholJacobian(pt.Lf[g], rowIdx, colIdx);
          Kc.submat(a, a, b, b) =
            cholCurvature(gradMatrix(gV, a, nT, rowIdx, colIdx), rowIdx, colIdx);
        }
        arma::mat Jf = J.submat(f, f);
        arma::vec gL = Jf.t()*gV.elem(f);
        // Negative Hessian of the likelihood with respect to the factors
        arma::mat negH = Jf.t()*AI.submat(f, f)*Jf - Kc.submat(f, f);
        negH = 0.5*(negH+negH.t());
        arma::vec hVal;
        arma::mat hVec;
        if(arma::eig_sym(hVal, hVec, negH)){
          // A direction of negative curvature has the sign of its eigenvalue
          // reversed, so the step is always one of ascent
          hVal = arma::abs(hVal);
          double floorVal = 1.0e-10*std::max(hVal.max(), 1.0e-300);
          hVal.elem(find(hVal<floorVal)).fill(floorVal);
          arma::vec step = hVec*((hVec.t()*gL)/hVal);
          double alpha = 1.0;
          for(int h=0; h<10; ++h){
            next.Lf = pt.Lf;
            arma::uword pos = 0;
            for(arma::uword g : freeMat){
              for(arma::uword q=0; q<nPar; ++q){
                next.Lf[g](rowIdx(q), colIdx(q)) += alpha*step(pos);
                ++pos;
              }
            }
            next.V.assign(nMat, arma::mat());
            for(arma::uword g=0; g<nMat; ++g){
              next.V[g] = next.Lf[g]*next.Lf[g].t();
            }
            if(evaluate(next) && (next.LL>=pt.LL)){
              accepted = true;
              break;
            }
            alpha *= 0.5;
          }
        }
      }
      if(!accepted){
        // The EM update V + V*S*V/r, with S the score matrix of V and r the
        // rank of its term, or an upper bound for it, which only shortens the
        // step and keeps V positive semidefinite. Ve is held at its floor. If
        // the result cannot be evaluated, the step is halved toward the
        // current point, which every point between the two keeps positive
        // semidefinite.
        std::vector<arma::mat> target(nMat);
        for(arma::uword g : freeMat){
          arma::mat S = gradMatrix(gV, g*nPar, nT, rowIdx, colIdx);
          S = S + S.t(); // the diagonal doubles, as gradParam halves it
          double r = (g<k) ? rankBound(g) : double(n);
          arma::mat Vn = pt.V[g] + pt.V[g]*S*pt.V[g]/r;
          target[g] = 0.5*(Vn+Vn.t());
        }
        bool emOk = false;
        double alpha = 1.0;
        for(int h=0; (h<=10) && !emOk; ++h){
          next.V = pt.V;
          next.Lf = pt.Lf;
          for(arma::uword g : freeMat){
            arma::mat Vg = pt.V[g] + alpha*(target[g]-pt.V[g]);
            if(g==k){
              Vg = floorVe(Vg);
            }
            next.Lf[g] = lowerChol(Vg);
            next.V[g] = next.Lf[g]*next.Lf[g].t();
          }
          emOk = evaluate(next);
          alpha *= 0.5;
        }
        if(!emOk){
          Rcpp::stop("%s: REML iteration left the parameter space", name);
        }
      }
      double change = next.LL - pt.LL;
      pt = next;
      if( (iter>nEM) && (std::abs(change)<tol) ){
        converged = true;
        break;
      }
    }
    if(!converged){
      iter = maxIter;
      Rf_warning("%s: did not converge in %d iterations", name, maxIter);
    }
    return iter;
  }

  // Generalized least squares estimates for the standardized traits, on the
  // orthonormal basis of the fixed effects
  arma::mat betaBasis(const MvMkPoint& pt) const{
    if(side==Side::record){
      arma::vec br = pt.Ainv*(pt.VinvX.t()*yv);
      return arma::reshape(br, p, nT);
    }
    return arma::reshape(pt.sol.head(p*nT), p, nT);
  }

  // BLUP of term i's effects for the standardized traits and the scaled
  // markers, U_i = M_i'*(P*y)*Vu_i on the record side and A_i*L_i' on the
  // coefficient side
  arma::mat blupScaled(const MvMkPoint& pt, arma::uword i) const{
    if(!haveCodes){
      Rcpp::stop("the BLUPs need the marker codes, which were not given");
    }
    if(held(i)){
      return arma::mat(Ms(i).n_cols, nT, arma::fill::zeros);
    }
    if(side==Side::record){
      return Ms(i).t()*pt.A*pt.V[i];
    }
    arma::uword len = segLen(i+1);
    if(len==0){
      return arma::mat(0, nT);
    }
    arma::uword off = nT*segStart(i+1);
    arma::mat S = arma::reshape(pt.sol.subvec(off, off+nT*len-1), len, nT);
    return S*pt.Lf[i].t();
  }

  // The mean diagonal each term was scaled by
  double termScale(arma::uword i) const{
    return held(i) ? 1.0 : md(i);
  }

private:
  Side side;
  arma::mat Y;   // the standardized responses with the least squares fit off
  arma::mat Xr;
  arma::field<arma::mat> Ms; // the scaled marker matrices
  bool haveCodes = true;     // false when built from the G_i alone
  arma::vec md, rankBound;
  arma::uvec held;
  arma::uvec rowIdx, colIdx;
  arma::uword nPar = 0;
  // Record side
  arma::field<arma::mat> G;
  arma::mat Xf;  // I (x) Xr
  arma::vec yv;  // vec(Y)
  // Coefficient side
  arma::uword nSeg = 0, nCoef = 0;
  arma::uvec segStart, segLen;
  arma::mat T, TtT, TtY;
  arma::field<arma::mat> GG; // T'*G_i*T

  // Sets up the fixed effects and the responses. A least squares fit is
  // taken off Y and added back to the fixed effects, and the traits are
  // standardized, for the reasons given for the single kernel solver; see
  // traitScale.
  void setResponses(const arma::mat& Yin, const arma::mat& X){
    n = Yin.n_rows;
    nT = Yin.n_cols;
    N = n*nT;
    fd = fixedDesign(X);
    if(fd.keep.n_elem==0){
      Rcpp::stop("X must have at least one nonzero column");
    }
    Xr = fd.Q;
    p = Xr.n_cols;
    if(p>=n){
      Rcpp::stop("the rank of X must be smaller than the number of records");
    }
    arma::mat Ys = removeLeastSquares(Xr, Yin, B0);
    ts = traitScale(Ys);
    Y = ts.Y;
    // The smallest eigenvalue Ve may take, as in the single trait solver: a
    // small fraction of the least residual variance of a standardized trait.
    // Without it a response that the markers fit almost exactly drives Ve, and
    // with it V, toward singular, where the likelihood cannot be evaluated.
    arma::rowvec vyFloor = arma::sum(arma::square(Y), 0)/double(n-p);
    double vyMin = vyFloor.min();
    if(!(vyMin>0.0) || !std::isfinite(vyMin)){
      vyMin = 1.0;
    }
    veFloor = 1e-8*vyMin;
  }

  // The factor of a block of the design in the trait dimension: the identity
  // for the fixed effects and L_i for term i
  arma::mat segFactor(const MvMkPoint& pt, arma::uword r) const{
    if(r==0){
      return arma::eye(nT, nT);
    }
    return pt.Lf[r-1];
  }

  // Unit matrix E_ab, symmetric
  arma::mat unitMatrix(arma::uword a, arma::uword b) const{
    arma::mat E(nT, nT, arma::fill::zeros);
    E(a,b) = 1.0;
    E(b,a) = 1.0;
    return E;
  }

  bool evaluateRecord(MvMkPoint& pt) const{
    pt.ok = false;
    arma::mat V(N, N, arma::fill::zeros);
    for(arma::uword a=0; a<nT; ++a){
      for(arma::uword b=0; b<nT; ++b){
        arma::mat blk(n, n, arma::fill::zeros);
        blk.diag().fill(pt.V[k](a,b));
        for(arma::uword i=0; i<k; ++i){
          if(held(i)) continue;
          blk += pt.V[i](a,b)*G(i);
        }
        V.submat(a*n, b*n, a*n+n-1, b*n+n-1) = blk;
      }
    }
    V = 0.5*(V+V.t()); // rounding only
    if(!arma::chol(pt.Lv, V, "lower")){
      return false;
    }
    V.reset();
    double logDetV = 2.0*arma::accu(arma::log(pt.Lv.diag()));
    const arma::solve_opts::opts fast = arma::solve_opts::fast;
    arma::mat Q = arma::solve(arma::trimatl(pt.Lv), Xf, fast);
    arma::vec z = arma::solve(arma::trimatl(pt.Lv), yv, fast);
    arma::mat AA = Q.t()*Q;
    AA = 0.5*(AA+AA.t()); // rounding only
    double logDetA;
    if(!invAndLogDet(pt.Ainv, logDetA, AA)){
      return false;
    }
    arma::vec r = z - Q*(pt.Ainv*(Q.t()*z));
    arma::vec Py = arma::solve(arma::trimatu(pt.Lv.t()), r, fast);
    pt.VinvX = arma::solve(arma::trimatu(pt.Lv.t()), Q, fast);
    pt.A = arma::reshape(Py, n, nT);
    double yPy = arma::dot(r, r);
    pt.LL = -0.5*(double(N-p*nT)*std::log(2.0*arma::datum::pi) + logDetV +
      logDetA + yPy);
    pt.ok = std::isfinite(pt.LL);
    return pt.ok;
  }

  void derivativesRecord(const MvMkPoint& pt, arma::vec& score,
                         arma::mat& AI) const{
    const arma::uword nTot = (k+1)*nPar;
    score.zeros(nTot);
    arma::mat P;
    {
      arma::mat Linv;
      if(!arma::inv(Linv, arma::trimatl(pt.Lv))){
        Rcpp::stop("could not invert the Cholesky factor of V");
      }
      P = Linv.t()*Linv;
    }
    P -= pt.VinvX*pt.Ainv*pt.VinvX.t();
    P = 0.5*(P+P.t()); // rounding only
    arma::mat Wv(N, nTot, arma::fill::zeros);
    for(arma::uword g=0; g<=k; ++g){
      if( (g<k) && held(g) ) continue;
      arma::mat GA = (g<k) ? arma::mat(G(g)*pt.A) : pt.A;
      arma::mat AGA = pt.A.t()*GA;
      for(arma::uword q=0; q<nPar; ++q){
        arma::uword a = rowIdx(q), b = colIdx(q);
        double quad = (a==b) ? AGA(a,a) : AGA(a,b)+AGA(b,a);
        double tr;
        if(g<k){
          tr = arma::accu(P.submat(a*n, b*n, a*n+n-1, b*n+n-1)%G(g));
        }else{
          tr = arma::trace(P.submat(a*n, b*n, a*n+n-1, b*n+n-1));
        }
        if(a!=b){
          tr *= 2.0;
        }
        score(g*nPar+q) = 0.5*(quad - tr);
        Wv.col(g*nPar+q) = arma::vectorise(GA*unitMatrix(a, b));
      }
    }
    AI = 0.5*(Wv.t()*(P*Wv));
    AI = 0.5*(AI+AI.t()); // rounding only
  }

  bool evaluateCoefficient(MvMkPoint& pt) const{
    pt.ok = false;
    arma::mat Le;
    if(!arma::chol(Le, pt.V[k], "lower")){
      return false;
    }
    arma::mat LeInv;
    if(!arma::inv(LeInv, arma::trimatl(Le))){
      return false;
    }
    pt.VeInv = LeInv.t()*LeInv;
    double logDetVe = 2.0*arma::accu(arma::log(Le.diag()));
    const arma::uword nC = nT*nCoef;
    arma::mat C(nC, nC);
    for(arma::uword r=0; r<nSeg; ++r){
      arma::mat Lr = segFactor(pt, r);
      for(arma::uword s=0; s<nSeg; ++s){
        if( (segLen(r)==0) || (segLen(s)==0) ) continue;
        arma::mat Ls = segFactor(pt, s);
        arma::mat blk = arma::kron(Lr.t()*pt.VeInv*Ls,
                                   TtT.submat(segStart(r), segStart(s),
                                              segStart(r)+segLen(r)-1,
                                              segStart(s)+segLen(s)-1));
        C.submat(nT*segStart(r), nT*segStart(s),
                 nT*segStart(r)+nT*segLen(r)-1,
                 nT*segStart(s)+nT*segLen(s)-1) = blk;
      }
    }
    for(arma::uword j=nT*p; j<nC; ++j){
      C(j,j) += 1.0;
    }
    C = 0.5*(C+C.t()); // rounding only
    if(!arma::chol(pt.Lc, C, "lower")){
      return false;
    }
    C.reset();
    double logDetC = 2.0*arma::accu(arma::log(pt.Lc.diag()));
    arma::vec rhs(nC);
    for(arma::uword r=0; r<nSeg; ++r){
      if(segLen(r)==0) continue;
      arma::mat blk = TtY.rows(segStart(r), segStart(r)+segLen(r)-1)*
        pt.VeInv*segFactor(pt, r);
      rhs.subvec(nT*segStart(r), nT*segStart(r)+nT*segLen(r)-1) =
        arma::vectorise(blk);
    }
    const arma::solve_opts::opts fast = arma::solve_opts::fast;
    pt.sol = arma::solve(arma::trimatu(pt.Lc.t()),
                         arma::solve(arma::trimatl(pt.Lc), rhs, fast), fast);
    arma::mat E = Y - multW(pt, pt.sol);
    pt.A = E*pt.VeInv;
    // y'*P*y = e'*R^-1*e + a'*a for the whitened effects a. This is equal to
    // y'*R^-1*e, but every term is a sum of squares, so nothing cancels: with
    // the residual variance close to zero, y'*R^-1*e would subtract numbers
    // amplified by Ve^-1 and lose the contribution of the random effects.
    double yPy = arma::accu(E%pt.A);
    if(nT*p<pt.sol.n_elem){
      yPy += arma::accu(arma::square(pt.sol.tail(pt.sol.n_elem-nT*p)));
    }
    pt.LL = -0.5*(double(N-p*nT)*std::log(2.0*arma::datum::pi) +
      double(n)*logDetVe + logDetC + yPy);
    pt.ok = std::isfinite(pt.LL);
    return pt.ok;
  }

  // W*s as an n by t matrix, sum_r T_r*S_r*L_r'
  arma::mat multW(const MvMkPoint& pt, const arma::vec& s) const{
    arma::mat out(n, nT, arma::fill::zeros);
    for(arma::uword r=0; r<nSeg; ++r){
      arma::uword len = segLen(r);
      if(len==0) continue;
      arma::uword off = nT*segStart(r);
      arma::mat S = arma::reshape(s.subvec(off, off+nT*len-1), len, nT);
      out += T.cols(segStart(r), segStart(r)+len-1)*S*segFactor(pt, r).t();
    }
    return out;
  }

  // W'*vec(Z) for an n by t matrix Z, segment by segment vec(T_r'*Z*L_r)
  arma::vec crossW(const MvMkPoint& pt, const arma::mat& Z) const{
    arma::vec out(nT*nCoef);
    for(arma::uword r=0; r<nSeg; ++r){
      arma::uword len = segLen(r);
      if(len==0) continue;
      arma::mat blk = T.cols(segStart(r), segStart(r)+len-1).t()*Z*
        segFactor(pt, r);
      out.subvec(nT*segStart(r), nT*segStart(r)+nT*len-1) =
        arma::vectorise(blk);
    }
    return out;
  }

  void derivativesCoefficient(const MvMkPoint& pt, arma::vec& score,
                              arma::mat& AI) const{
    const arma::uword nTot = (k+1)*nPar;
    score.zeros(nTot);
    arma::mat Cinv;
    {
      arma::mat Linv;
      if(!arma::inv(Linv, arma::trimatl(pt.Lc))){
        Rcpp::stop("could not invert the Cholesky factor of C");
      }
      Cinv = Linv.t()*Linv;
    }
    arma::mat Wv(N, nTot, arma::fill::zeros);
    for(arma::uword g=0; g<=k; ++g){
      if( (g<k) && held(g) ) continue;
      const arma::mat& GGg = (g<k) ? GG(g) : TtT;
      // Gam(r,s)(t,u) is the trace of the (t,u) block of the (r,s) block of
      // C^-1 against T_r'*G*T_s, which every element of the matrix needs
      arma::field<arma::mat> Gam(nSeg, nSeg);
      for(arma::uword r=0; r<nSeg; ++r){
        for(arma::uword s=0; s<nSeg; ++s){
          Gam(r,s).zeros(nT, nT);
          if( (segLen(r)==0) || (segLen(s)==0) ) continue;
          arma::mat GGrs = GGg.submat(segStart(r), segStart(s),
                                      segStart(r)+segLen(r)-1,
                                      segStart(s)+segLen(s)-1);
          for(arma::uword t=0; t<nT; ++t){
            arma::uword r0 = nT*segStart(r)+t*segLen(r);
            for(arma::uword u=0; u<nT; ++u){
              arma::uword s0 = nT*segStart(s)+u*segLen(s);
              Gam(r,s)(t,u) = arma::accu(Cinv.submat(r0, s0,
                                                     r0+segLen(r)-1,
                                                     s0+segLen(s)-1)%GGrs);
            }
          }
        }
      }
      arma::mat GA, AGA;
      if(g<k){
        arma::mat TA = Ms(g).t()*pt.A;
        AGA = TA.t()*TA;
        GA = Ms(g)*TA;
      }else{
        AGA = pt.A.t()*pt.A;
        GA = pt.A;
      }
      for(arma::uword q=0; q<nPar; ++q){
        arma::uword a = rowIdx(q), b = colIdx(q);
        arma::mat E = unitMatrix(a, b);
        double quad = (a==b) ? AGA(a,a) : AGA(a,b)+AGA(b,a);
        arma::mat F = pt.VeInv*E*pt.VeInv;
        double trWCW = 0.0;
        for(arma::uword r=0; r<nSeg; ++r){
          arma::mat Lr = segFactor(pt, r);
          for(arma::uword s=0; s<nSeg; ++s){
            trWCW += arma::accu((Lr.t()*F*segFactor(pt, s))%Gam(r,s));
          }
        }
        // tr(G) is n for every scaled term and for the identity
        double tr = arma::trace(pt.VeInv*E)*double(n) - trWCW;
        score(g*nPar+q) = 0.5*(quad - tr);
        Wv.col(g*nPar+q) = arma::vectorise(GA*E);
      }
    }
    // P*w = R^-1*(w - W*C^-1*W'*R^-1*w), column by column
    arma::mat PW(N, nTot, arma::fill::zeros);
    for(arma::uword j=0; j<nTot; ++j){
      if(!arma::any(Wv.col(j))) continue;
      arma::mat Wm = arma::reshape(Wv.col(j), n, nT);
      arma::mat RW = Wm*pt.VeInv;
      arma::vec x = Cinv*crossW(pt, RW);
      PW.col(j) = arma::vectorise((Wm - multW(pt, x))*pt.VeInv);
    }
    AI = 0.5*(Wv.t()*PW);
    AI = 0.5*(AI+AI.t()); // rounding only
  }
};

// Fits the multiple trait RR-BLUP form of the multikernel model, on the side
// requested as in fitMultiKernel. The inputs must already have been checked.
Rcpp::List fitMultiKernelMV(const arma::mat& Y, const arma::mat& X,
                            const arma::field<arma::mat>& Mlist,
                            int maxIter, double tol, int sideRequest,
                            const char* name, bool reportSide){
  const arma::uword n = Y.n_rows;
  const arma::uword nT = Y.n_cols;
  const arma::uword k = Mlist.n_elem;
  Side side;
  if(sideRequest==1){
    side = Side::record;
  }else if(sideRequest==2){
    side = Side::coefficient;
  }else{
    // The costs of useCoefficientSide with the records, the coefficients and
    // the parameters of the multiple trait model
    arma::uword c = X.n_cols;
    for(arma::uword i=0; i<k; ++i){
      c += Mlist(i).n_cols;
    }
    side = useCoefficientSide(n*nT, c*nT, (k+1)*nT*(nT+1)/2) ?
      Side::coefficient : Side::record;
  }
  MvMultiKernelREML reml(Y, X, side, Mlist);
  MvMkPoint pt;
  int iter = reml.fit(pt, maxIter, tol, name);

  // Back to the units of the traits and of the markers. With D = diag(sd),
  // the covariance matrices are D*V*D, divided for a term by the mean
  // diagonal it was scaled by, and the effects are scaled the same way.
  const arma::mat D = arma::diagmat(reml.ts.sd);
  Rcpp::List Vu(k), u(k);
  for(arma::uword i=0; i<k; ++i){
    double s = reml.termScale(i);
    Vu[i] = Rcpp::wrap(arma::mat(D*pt.V[i]*D/s));
    u[i] = Rcpp::wrap(arma::mat(reml.blupScaled(pt, i)*D/std::sqrt(s)));
  }
  arma::mat Ve = D*pt.V[k]*D;
  arma::mat beta = fixedEffectsFromBasis(reml.fd, reml.betaBasis(pt)*D+reml.B0,
                                         X.n_cols);
  double p = double(reml.fd.Q.n_cols);
  double LL = pt.LL - (double(n)-p)*arma::accu(arma::log(reml.ts.sd)) -
    double(nT)*reml.fd.logDetR;
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                                      Rcpp::Named("Ve")=Ve,
                                      Rcpp::Named("beta")=beta,
                                      Rcpp::Named("u")=u,
                                      Rcpp::Named("LL")=LL,
                                      Rcpp::Named("iter")=iter);
  if(reportSide){
    out.push_back((side==Side::record) ? 1 : 2, "side");
  }
  return out;
}

// The checks of checkResponse for a response with any number of columns
void checkResponseMV(const arma::mat& Y, const arma::mat& X, const char* name){
  if(X.n_rows!=Y.n_rows){
    Rcpp::stop("%s: X and y must have the same number of rows", name);
  }
  if(!Y.is_finite() || !X.is_finite()){
    Rcpp::stop("%s: y and X must not contain missing or infinite values",
               name);
  }
}

} // namespace

// Variance components of the RR-BLUP form of the multikernel model, for one
// trait or several, from each term's G_i = M_i*M_i' rather than from M_i.
// The fast mode of the multikernel models forms the G_i of a subset of
// records itself when there are at least as many markers as records, and the
// record side works with G_i alone, so passing them on saves finding a thin
// factor of each, which costs a decomposition of every G_i with its
// eigenvectors. Vu and Ve are returned as solveRRBLUPMK returns them, without
// the BLUPs. rank is an upper bound for the rank of each G_i, and the G_i
// are moved out of Glist. The inputs must already have been checked.
Rcpp::List remlMultiKernelGram(const arma::mat& Y, const arma::mat& X,
                               arma::field<arma::mat>& Glist,
                               const arma::vec& rank, int maxIter,
                               const char* name){
  const arma::uword k = Glist.n_elem;
  if(Y.n_cols>1){
    MvMultiKernelREML reml(Y, X, Glist, rank);
    MvMkPoint pt;
    int iter = reml.fit(pt, maxIter, rrblupMKTol, name);
    // Back to the units of the traits and of the markers, as in
    // fitMultiKernelMV
    const arma::mat D = arma::diagmat(reml.ts.sd);
    Rcpp::List Vu(k);
    for(arma::uword i=0; i<k; ++i){
      Vu[i] = Rcpp::wrap(arma::mat(D*pt.V[i]*D/reml.termScale(i)));
    }
    arma::mat Ve = D*pt.V[k]*D;
    return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                              Rcpp::Named("Ve")=Ve,
                              Rcpp::Named("iter")=iter);
  }
  arma::vec yv = Y.col(0);
  MultiKernelREML reml(yv, X, Glist, rank);
  RemlPoint pt;
  int iter = reml.fit(pt, maxIter, rrblupMKTol, name);
  arma::mat Vu = pt.theta.head(k);
  arma::mat Ve(1, 1);
  Ve(0,0) = pt.theta(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("iter")=iter);
}

// Fits either form of the multikernel model on a chosen side, for testing.
// side is 0 to choose automatically, 1 for the record side and 2 for the
// coefficient side, and the side used is returned as an extra element.
// [[Rcpp::export(.solveMultiKernel)]]
Rcpp::List solveMultiKernel(const arma::mat& y, const arma::mat& X,
                            const arma::field<arma::mat>& Zlist,
                            const arma::field<arma::mat>& Klist, bool rrForm,
                            int maxIter, double tol, int side){
  const char* name = rrForm ? "solveRRBLUPMK" : "solveMKM";
  if(rrForm && (y.n_cols>1)){
    checkResponseMV(y, X, name);
    checkMlist(y, Zlist, name);
    return fitMultiKernelMV(y, X, Zlist, maxIter, tol, side, name, true);
  }
  checkResponse(y, X, name);
  if(rrForm){
    checkMlist(y, Zlist, name);
  }else{
    checkZKlist(y, Zlist, Klist, name);
  }
  return fitMultiKernel(y, X, Zlist, Klist, rrForm, maxIter, tol, side,
                        name, true);
}

//' @title Solve Multikernel RR-BLUP
//'
//' @description
//' Solves a mixed model with multiple random effects, where random term i
//' contributes \code{Vu[i]*Mlist[[i]]\%*\%t(Mlist[[i]])} to the variance of
//' the records and its effects are the marker effects of \code{Mlist[[i]]}.
//' The variance components are estimated by average information REML
//' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}, and the
//' matrices in \code{Mlist} may have different numbers of columns. When
//' \code{y} has more than one column, a multiple trait model is fitted, in
//' which each term's effects and the errors have unstructured covariance
//' matrices among the traits.
//'
//' @details
//' The model is fitted in whichever of two equivalent ways is cheaper. With
//' fewer records than marker effects it works with the variance matrix of the
//' records, which has one row per record. With fewer marker effects than
//' records it works with Henderson's mixed model equations, which have one
//' row per effect. Both give the same estimates; only the time and memory
//' they need differ.
//'
//' Iteration starts with half of the residual variance from a least
//' squares fit of X assigned to the error and half shared equally among
//' the random terms. A step that would lower the likelihood is halved, and
//' an EM-REML step is used if halving does not help. Variance components
//' are kept non-negative: a component that reaches zero stays there for as
//' long as the likelihood would fall by increasing it. Iteration stops when
//' the change in log-likelihood and the relative change in the variance
//' components are both below 1e-6, and a warning is given if
//' \code{maxIter} is reached first.
//'
//' The multiple trait model needs every trait recorded on every record. Its
//' covariance matrices start with half of each trait's residual variance from
//' the least squares fit assigned to the error and half shared equally among
//' the random terms, with the correlation among the least squares residuals
//' for every matrix. Five EM iterations follow, and then Newton steps in the
//' Cholesky factors of the covariance matrices, as described for
//' \code{\link{solveRRBLUPMV}}, which keeps every matrix positive
//' semidefinite. Iteration stops when the restricted log-likelihood changes
//' by less than 1e-6.
//'
//' @param y a matrix with n rows and 1 column, or one column per trait for a
//' multiple trait model
//' @param X a matrix with n rows and x columns. Columns that are linearly
//' dependent on the others are dropped and get a fixed effect of zero. An
//' intercept, a column whose values are all the same, is always kept.
//' @param Mlist a list of M matrices
//' @param maxIter maximum number of iterations
//'
//' @return A list with the variance components \code{Vu} (one per matrix in
//' \code{Mlist}) and \code{Ve}, the fixed effects \code{beta}, a list
//' \code{u} of marker effects for each matrix in \code{Mlist}, the
//' restricted log-likelihood \code{LL} at the estimates, and the number of
//' iterations \code{iter}. \code{LL} is
//' \eqn{-\frac{1}{2}[(n-p)\log(2\pi) + \log|V| + \log|X'V^{-1}X| + y'Py]},
//' where p is the rank of X; it excludes the \eqn{\log|X'X|} term.
//'
//' For a multiple trait model with q traits, \code{Vu} is instead a list of
//' q by q covariance matrices, one per matrix in \code{Mlist}, \code{Ve} is a
//' q by q matrix, \code{beta} has one column per trait, and each element of
//' \code{u} has one column per trait. \code{LL} is then that of
//' \eqn{y=vec(Y)}, with nq records and p equal to q times the rank of X.
//'
//' @references
//' \insertAllCited{}
//'
//' @examples
//' set.seed(1)
//' n = 50
//' M1 = matrix(rnorm(n*20), n)
//' M2 = matrix(rnorm(n*10), n)
//' y = M1%*%rnorm(20, sd=0.3) + rnorm(n)
//' ans = solveRRBLUPMK(y, matrix(1, n, 1), list(M1, M2))
//' ans$Vu
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMK(const arma::mat& y, const arma::mat& X,
                         const arma::field<arma::mat>& Mlist,
                         int maxIter=40){
  const char* name = "solveRRBLUPMK";
  if(y.n_cols>1){
    checkResponseMV(y, X, name);
    checkMlist(y, Mlist, name);
    return fitMultiKernelMV(y, X, Mlist, maxIter, rrblupMKTol, 0, name,
                            false);
  }
  checkResponse(y, X, name);
  checkMlist(y, Mlist, name);
  arma::field<arma::mat> none;
  return fitMultiKernel(y, X, Mlist, none, true, maxIter, rrblupMKTol, 0,
                        name, false);
}

//' @title Solve Multikernel Model
//'
//' @description
//' Solves a univariate mixed model with multiple random effects,
//' \eqn{y = X\beta + \sum_i Z_i u_i + e}, where \eqn{u_i} has covariance
//' \eqn{K_i V_{u_i}} and \eqn{e} has covariance \eqn{I V_e}. The variance
//' components are estimated by average information REML
//' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}.
//'
//' @details
//' The model is fitted in whichever of two equivalent ways is cheaper. When
//' the random terms have more levels in total than there are records, it
//' works with the variance matrix of the records,
//' \eqn{V = \sum_i V_{u_i} Z_i K_i Z_i' + I V_e}. Otherwise it works with
//' Henderson's mixed model equations, writing each \eqn{K_i} as
//' \eqn{F_i F_i'} over its positive eigenvalues, so a singular \eqn{K_i} is
//' allowed. Both give the same estimates; only the time and memory they need
//' differ. Each \eqn{K_i} must be positive semidefinite.
//'
//' Iteration starts with half of the residual variance from a least
//' squares fit of X assigned to the error and half shared equally among
//' the random terms. A step that would lower the likelihood is halved, and
//' an EM-REML step is used if halving does not help. Variance components
//' are kept non-negative: a component that reaches zero stays there for as
//' long as the likelihood would fall by increasing it. Iteration stops when
//' the change in log-likelihood and the relative change in the variance
//' components are both below \code{tol}, and a warning is given if
//' \code{maxIter} is reached first.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns. Columns that are linearly
//' dependent on the others are dropped and get a fixed effect of zero. An
//' intercept, a column whose values are all the same, is always kept.
//' @param Zlist a list of Z matrices
//' @param Klist a list of K matrices
//' @param maxIter maximum number of iterations
//' @param tol tolerance for convergence
//'
//' @return A list with the variance components \code{Vu} (one per random
//' term) and \code{Ve}, the fixed effects \code{beta}, a list \code{u} of
//' BLUPs for each random term, the restricted log-likelihood \code{LL} at
//' the estimates, and the number of iterations \code{iter}. \code{LL} is
//' \eqn{-\frac{1}{2}[(n-p)\log(2\pi) + \log|V| + \log|X'V^{-1}X| + y'Py]},
//' where p is the rank of X; it excludes the \eqn{\log|X'X|} term.
//'
//' @references
//' \insertAllCited{}
//'
//' @examples
//' set.seed(1)
//' n = 50
//' M = matrix(rnorm(n*20), n)
//' K = tcrossprod(M)/20
//' y = M%*%rnorm(20, sd=0.3) + rnorm(n)
//' Z = diag(n)
//' ans = solveMKM(y, matrix(1, n, 1), list(Z, Z), list(K, diag(n)))
//' ans$Vu
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMKM(const arma::mat& y, const arma::mat& X,
                    const arma::field<arma::mat>& Zlist,
                    const arma::field<arma::mat>& Klist,
                    int maxIter=40, double tol=1e-4){
  const char* name = "solveMKM";
  checkResponse(y, X, name);
  checkZKlist(y, Zlist, Klist, name);
  return fitMultiKernel(y, X, Zlist, Klist, false, maxIter, tol, 0,
                        name, false);
}
