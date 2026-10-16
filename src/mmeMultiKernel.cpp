#include "mme.h"

// Single trait mixed models with several random terms,
//
//   y = X*beta + sum_i Z_i*u_i + e,  u_i ~ N(0, K_i*Vu_i),  e ~ N(0, I*Ve)
//
// fitted by average information REML (Gilmour, Thompson and Cullis 1995;
// Johnson and Thompson 1995). Everything is done on the marginal variance
//
//   V = sum_i Vu_i*G_i + Ve*I,  G_i = Z_i*K_i*Z_i'
//
// rather than on Henderson's mixed model equations, because the number of
// records is usually far smaller than the number of markers in the RR-BLUP
// form, and V has the dimension of the records whatever the kernels are.
//
// With P = V^-1 - V^-1*X*(X'*V^-1*X)^-1*X'*V^-1, the derivative of the REML
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
// Each G_i is held as a product L_i*R_i' of two thin matrices, L_i = Z_i*K_i
// and R_i = Z_i, or L_i = R_i = M_i in the RR-BLUP form. When those have
// fewer columns than there are records the n by n matrix G_i is never
// formed: G_i*v is L_i*(R_i'*v) and tr(P*G_i) = tr((P*L_i)*R_i') is the sum
// of the elementwise product of P*L_i and R_i. Otherwise forming G_i once is
// cheaper than repeating those products in every iteration, and it is kept.

namespace {

// Fixed tolerance for solveRRBLUPMK, which takes no tol argument. It is
// applied the same way as tol in solveMKM.
const double rrblupMKTol = 1e-6;

// Most times the step length can be halved before falling back to an EM step
const int maxHalving = 10;

// Everything the iterations need at one value of the variance components
struct RemlPoint {
  arma::vec theta; // Vu_1..Vu_k, then Ve
  arma::mat P;
  arma::vec Py;
  arma::mat VinvX;
  arma::mat Ainv;  // (X'*V^-1*X)^-1
  double LL = 0.0;
  bool ok = false;
};

class MultiKernelREML {
public:
  MultiKernelREML(const arma::vec& y_, const arma::mat& X,
                  const std::vector<const arma::mat*>& L_,
                  const std::vector<const arma::mat*>& R_)
    : y(y_), L(L_), R(R_){
    n = y.n_elem;
    k = L.size();
    xKeep = independentColumns(X);
    if(xKeep.n_elem==0){
      Rcpp::stop("X must have at least one nonzero column");
    }
    Xr = X.cols(xKeep);
    p = Xr.n_cols;
    if(p>=n){
      Rcpp::stop("the rank of X must be smaller than the number of records");
    }
    G.set_size(k);
    dense.assign(k, false);
    meanDiag.set_size(k);
    rankBound.set_size(k);
    for(arma::uword i=0; i<k; ++i){
      const arma::mat& Li = *L[i];
      const arma::mat& Ri = *R[i];
      meanDiag(i) = arma::accu(Li%Ri)/double(n);
      rankBound(i) = double(std::min(n, Li.n_cols));
      if(Li.n_cols>=n){
        dense[i] = true;
        G(i) = Li*Ri.t();
        G(i) = 0.5*(G(i)+G(i).t()); // rounding only
      }
    }
  }

  // Number of records, terms and the rank of X
  arma::uword n, k, p;
  arma::uvec xKeep;

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

  // Factorizes V at pt.theta and fills in P and the REML log-likelihood
  //
  //   l = -0.5*[(n-p)*log(2*pi) + log|V| + log|X'*V^-1*X| + y'*P*y]
  //
  // Returns false if V or X'*V^-1*X is not positive definite.
  bool evaluate(RemlPoint& pt) const{
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

  // Score vector and average information matrix at an evaluated point
  void derivatives(const RemlPoint& pt, arma::vec& score,
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
      if(!accepted){
        // No step improves the likelihood, so this is a stationary point to
        // within the precision of the arithmetic
        converged = true;
        break;
      }
      double dLL = trial.LL - pt.LL;
      double dTheta = relativeChange(trial.theta, pt.theta);
      pt = trial;
      if(std::abs(dLL)<tol && dTheta<tol){
        converged = true;
        break;
      }
    }
    if(!converged){
      Rf_warning("%s: did not converge in %d iterations", name, maxIter);
    }
    return iter;
  }

  // Generalized least squares estimates, with zeros for any columns of X
  // that were dropped for being linearly dependent on earlier ones
  arma::mat beta(const RemlPoint& pt, arma::uword nColX) const{
    arma::vec br = pt.Ainv*(pt.VinvX.t()*y);
    arma::mat b(nColX, 1, arma::fill::zeros);
    b.elem(xKeep) = br;
    return b;
  }

  // BLUP of the term's effects, u_i = Vu_i*K_i*Z_i'*P*y = Vu_i*L_i'*P*y
  arma::mat blup(const RemlPoint& pt, arma::uword i) const{
    arma::mat u = pt.theta(i)*(L[i]->t()*pt.Py);
    return u;
  }

private:
  const arma::vec& y;
  arma::mat Xr;
  std::vector<const arma::mat*> L, R;
  arma::field<arma::mat> G;
  std::vector<bool> dense;
  arma::vec meanDiag, rankBound;
  arma::uvec estimable;
  double vy = 1.0, veFloor = 0.0;
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

Rcpp::List packResult(const MultiKernelREML& reml, const RemlPoint& pt,
                      arma::uword nColX, int iter){
  arma::uword k = reml.k;
  arma::field<arma::mat> u(k);
  for(arma::uword i=0; i<k; ++i){
    u(i) = reml.blup(pt, i);
  }
  arma::mat Vu = pt.theta.head(k);
  arma::mat Ve(1, 1);
  Ve(0,0) = pt.theta(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=reml.beta(pt, nColX),
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=pt.LL,
                            Rcpp::Named("iter")=iter);
}

} // namespace

//' @title Solve Multikernel RR-BLUP
//'
//' @description
//' Solves a univariate mixed model with multiple random effects, where
//' random term i contributes \code{Vu[i]*Mlist[[i]]\%*\%t(Mlist[[i]])} to the
//' variance of the records and its effects are the marker effects of
//' \code{Mlist[[i]]}. The variance components are estimated by average
//' information REML \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}
//' working on the variance matrix of the records, so the matrices in
//' \code{Mlist} may have different numbers of columns and there may be more
//' markers than records.
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
  arma::uword k = Mlist.n_elem;
  if(k==0){
    Rcpp::stop("%s: Mlist must contain at least one matrix", name);
  }
  std::vector<const arma::mat*> L(k);
  for(arma::uword i=0; i<k; ++i){
    if(Mlist(i).n_rows!=y.n_rows){
      Rcpp::stop("%s: every matrix in Mlist must have one row per record",
                 name);
    }
    if(!Mlist(i).is_finite()){
      Rcpp::stop("%s: Mlist must not contain missing or infinite values",
                 name);
    }
    L[i] = &Mlist(i);
  }
  arma::vec yv = y.col(0);
  MultiKernelREML reml(yv, X, L, L);
  RemlPoint pt;
  int iter = reml.fit(pt, maxIter, rrblupMKTol, name);
  return packResult(reml, pt, X.n_cols, iter);
}

//' @title Solve Multikernel Model
//'
//' @description
//' Solves a univariate mixed model with multiple random effects,
//' \eqn{y = X\beta + \sum_i Z_i u_i + e}, where \eqn{u_i} has covariance
//' \eqn{K_i V_{u_i}} and \eqn{e} has covariance \eqn{I V_e}. The variance
//' components are estimated by average information REML
//' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR} working
//' on the variance matrix of the records,
//' \eqn{V = \sum_i V_{u_i} Z_i K_i Z_i' + I V_e}.
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
  arma::uword k = Zlist.n_elem;
  if(k==0){
    Rcpp::stop("%s: Zlist must contain at least one matrix", name);
  }
  if(Klist.n_elem!=k){
    Rcpp::stop("%s: Zlist and Klist must have the same length", name);
  }
  // Z_i*K_i is formed once. K_i is symmetrized first so that Z_i*K_i*Z_i'
  // is symmetric and the BLUP K_i*Z_i'*P*y can be read off as (Z_i*K_i)'*P*y.
  arma::field<arma::mat> ZK(k);
  std::vector<const arma::mat*> L(k), R(k);
  for(arma::uword i=0; i<k; ++i){
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
    ZK(i) = Z*(0.5*(K+K.t()));
    L[i] = &ZK(i);
    R[i] = &Zlist(i);
  }
  arma::vec yv = y.col(0);
  MultiKernelREML reml(yv, X, L, R);
  RemlPoint pt;
  int iter = reml.fit(pt, maxIter, tol, name);
  return packResult(reml, pt, X.n_cols, iter);
}
