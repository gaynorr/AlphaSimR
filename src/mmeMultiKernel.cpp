#include "mme.h"

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
  arma::mat P;
  arma::vec Py;
  arma::mat VinvX;
  arma::mat Ainv;  // (X'*V^-1*X)^-1
  // Coefficient side. A term whose variance is zero has no effects, so it is
  // left out of the equations rather than given an infinite precision.
  std::vector<Segment> segs; // the fixed effects, then each term in the equations
  std::vector<long> termSeg; // segment of each term, or -1 if left out
  arma::uvec idx;            // columns of T in the equations
  arma::mat Cinv;
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
    n = y.n_elem;
    k = L.size();
    // Columns of X that depend on the others are dropped, and the fit uses
    // Xr, an orthonormal basis for the span of the rest, so a badly scaled X
    // costs no precision; see fixedDesign
    fd = fixedDesign(X);
    if(fd.keep.n_elem==0){
      Rcpp::stop("X must have at least one nonzero column");
    }
    Xr = fd.Q;
    p = Xr.n_cols;
    if(p>=n){
      Rcpp::stop("the rank of X must be smaller than the number of records");
    }
    // A least squares fit is taken off y and added back to the fixed effects
    // by beta(), which keeps a large mean from swamping the sums of squares
    arma::mat B0;
    y = removeLeastSquares(Xr, y, B0).col(0);
    b0 = B0.col(0);
    meanDiag.set_size(k);
    rankBound.set_size(k);
    for(arma::uword i=0; i<k; ++i){
      const arma::mat& Li = *L[i];
      const arma::mat& Ri = *R[i];
      meanDiag(i) = arma::accu(Li%Ri)/double(n);
      rankBound(i) = double(std::min(n, Li.n_cols));
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
  // diagonal explains nothing and is held at zero throughout.
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
      estimable(i) = (meanDiag(i)>0.0) && std::isfinite(meanDiag(i));
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
  // together with what the derivatives need. Returns false if V, or the
  // matrix that stands in for it, is not positive definite.
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

  // First-order optimality on the boundary. A component held at zero is
  // optimal only if moving it up would not raise the likelihood, which to
  // first order is a gain of score*Vu for a move to Vu. The gain from a move
  // to a typical variance must be below a small multiple of the tolerance.
  bool boundaryOptimal(const RemlPoint& pt, const arma::vec& score,
                       double tol) const{
    for(arma::uword i=0; i<k; ++i){
      if(estimable(i) && (pt.theta(i)==0.0) &&
         (score(i)*typicalVu(i) > 10.0*tol)){
        return false;
      }
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
        // A small change is only convergence if no component held at zero
        // could raise the likelihood by moving off it
        derivatives(pt, score, AI);
        if(boundaryOptimal(pt, score, tol)){
          converged = true;
          break;
        }
      }else if(boundaryOptimal(pt, score, tol)){
        // No step improves the likelihood and none would from the boundary,
        // so this is a stationary point to within the precision of the
        // arithmetic
        converged = true;
        break;
      }
      // A component at zero has a positive score but neither step has moved
      // it, and the EM step cannot, because it scales the change by the
      // square of the component. A projected gradient step can.
      RemlPoint released;
      if(!releaseStep(pt, score, released)){
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
  arma::uvec estimable;
  double vy = 1.0, veFloor = 0.0;
  // Record side
  arma::field<arma::mat> G;
  std::vector<bool> dense;
  // Coefficient side
  arma::mat T, TtT;
  arma::vec Tty;
  arma::uvec termStart, termLen;

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
    arma::mat Vinv;
    double logDetV, logDetA;
    if(!invAndLogDet(Vinv, logDetV, V)){
      return false;
    }
    V.reset();
    pt.VinvX = Vinv*Xr;
    arma::mat A = Xr.t()*pt.VinvX;
    A = 0.5*(A+A.t()); // rounding only
    if(!invAndLogDet(pt.Ainv, logDetA, A)){
      return false;
    }
    pt.P = Vinv - pt.VinvX*pt.Ainv*pt.VinvX.t();
    pt.P = 0.5*(pt.P+pt.P.t()); // rounding only
    pt.Py = pt.P*y;
    double yPy = arma::dot(y, pt.Py);
    pt.LL = -0.5*(double(n-p)*std::log(2.0*arma::datum::pi) + logDetV +
      logDetA + yPy);
    pt.ok = std::isfinite(pt.LL);
    return pt.ok;
  }

  void derivativesRecord(const RemlPoint& pt, arma::vec& score,
                         arma::mat& AI) const{
    score.set_size(k+1);
    arma::mat W(n, k+1);
    for(arma::uword i=0; i<k; ++i){
      double trPG;
      if(dense[i]){
        W.col(i) = G(i)*pt.Py;
        trPG = arma::accu(pt.P%G(i));
      }else{
        W.col(i) = (*L[i])*(R[i]->t()*pt.Py);
        trPG = arma::accu((pt.P*(*L[i]))%(*R[i]));
      }
      score(i) = 0.5*(arma::dot(pt.Py, W.col(i)) - trPG);
    }
    W.col(k) = pt.Py;
    score(k) = 0.5*(arma::dot(pt.Py, pt.Py) - arma::trace(pt.P));
    AI = 0.5*(W.t()*(pt.P*W));
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
    double logDetC;
    if(!invAndLogDet(pt.Cinv, logDetC, C)){
      return false;
    }
    C.reset();
    arma::vec rhs = Tty.elem(pt.idx)/ve;
    pt.sol = pt.Cinv*rhs;
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
        double trCii = arma::trace(pt.Cinv.submat(s.pos, s.pos,
                                                  s.pos+s.len-1,
                                                  s.pos+s.len-1));
        trPG = double(s.len)/th - trCii/(th*th);
      }else{
        // A term left out of the equations: tr(P*G_i) is
        // tr(W_i'*W_i)/Ve - tr(B*C^-1*B')/Ve^2 for B = W_i'*T_a
        arma::uvec rows = arma::regspace<arma::uvec>(a, b);
        arma::mat B = TtT.submat(rows, pt.idx);
        trPG = arma::trace(TtT.submat(a, a, b, b))/ve -
          arma::accu((B*pt.Cinv)%B)/(ve*ve);
      }
      score(i) = 0.5*(quad - trPG);
    }
    W.col(k) = pt.e/ve;
    double trP = double(n)/ve -
      arma::accu(TtT.submat(pt.idx, pt.idx)%pt.Cinv)/(ve*ve);
    score(k) = 0.5*(arma::dot(pt.e, pt.e)/(ve*ve) - trP);
    // P*W = (W - T_a*C^-1*T_a'*W/Ve)/Ve
    arma::mat PW = (W - multT(pt, pt.Cinv*crossT(pt, W))/ve)/ve;
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

} // namespace

// Fits either form of the multikernel model on a chosen side, for testing.
// side is 0 to choose automatically, 1 for the record side and 2 for the
// coefficient side, and the side used is returned as an extra element.
// [[Rcpp::export(.solveMultiKernel)]]
Rcpp::List solveMultiKernel(arma::mat& y, arma::mat& X,
                            arma::field<arma::mat>& Zlist,
                            arma::field<arma::mat>& Klist, bool rrForm,
                            int maxIter, double tol, int side){
  const char* name = rrForm ? "solveRRBLUPMK" : "solveMKM";
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
//' Solves a univariate mixed model with multiple random effects, where
//' random term i contributes \code{Vu[i]*Mlist[[i]]\%*\%t(Mlist[[i]])} to the
//' variance of the records and its effects are the marker effects of
//' \code{Mlist[[i]]}. The variance components are estimated by average
//' information REML \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR},
//' and the matrices in \code{Mlist} may have different numbers of columns.
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
//' @param y a matrix with n rows and 1 column
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
Rcpp::List solveRRBLUPMK(arma::mat& y, arma::mat& X,
                         arma::field<arma::mat>& Mlist,
                         int maxIter=40){
  const char* name = "solveRRBLUPMK";
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
Rcpp::List solveMKM(arma::mat& y, arma::mat& X,
                    arma::field<arma::mat>& Zlist,
                    arma::field<arma::mat>& Klist,
                    int maxIter=40, double tol=1e-4){
  const char* name = "solveMKM";
  checkResponse(y, X, name);
  checkZKlist(y, Zlist, Klist, name);
  return fitMultiKernel(y, X, Zlist, Klist, false, maxIter, tol, 0,
                        name, false);
}
