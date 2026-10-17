#include "alphasimr.h"

#if !defined(ARMA_BLAS_CAPITALS)
#define arma_dsyevr dsyevr
#else
#define arma_dsyevr DSYEVR
#endif

extern "C"
void arma_fortran(arma_dsyevr)(char* JOBZ, char* RANGE, char* UPLO, long long int* N, double* A, long long int* LDA, double* VL,
                  double* VU, long long int* IL, long long int* IU, double* ABSTOL, long long int* M, double* W, double* Z,
                  long long int* LDZ, long long int* ISUPPZ, double* WORK, long long int* LWORK, long long int* IWORK,
                  long long int* LIWORK, long long int* INFO);

#include "mme.h"
#include "packedGv.h"

// Solvers defined in mmeMultiTrait.cpp and mmeMultiKernel.cpp. The defaults
// repeat those of the definitions, which Rcpp reads to build the R
// interface, so they are declared here rather than in mme.h, where the
// defining files would see them twice.
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X,
                         const arma::mat& M, int maxIter=1000,
                         double tol=1e-6);
Rcpp::List solveRRBLUPMK(const arma::mat& y, const arma::mat& X,
                         const arma::field<arma::mat>& Mlist,
                         int maxIter=40);

// Factorization of a mixed model coefficient matrix.
//
// The EM algorithm needs only two things from the coefficient matrix C: the
// solution of C*x=b, and the trace of the diagonal block of C^-1 belonging to
// each random effect. Neither requires the full inverse.
//
// C is symmetric and, with a positive value on the diagonal of every random
// effect block, positive definite, so C=L*L.t(). Writing Linv for inv(L),
//
//   C^-1     = Linv.t()*Linv
//   C^-1*b   = Linv.t()*(Linv*b)
//   diag(C^-1)_i = sum of squares of column i of Linv
//
// so the trace of a diagonal block is a sum of squares over the matching
// columns. Building Linv costs a Cholesky factorization plus a triangular
// inverse, about two thirds of what forming the full inverse costs, and both
// quantities above are then quadratic rather than cubic in the matrix size.
//
// A matrix that is not positive definite, which happens when X is rank
// deficient, falls back to the general inverse so that behavior is unchanged.
//
// Only the EM algorithm asks for a trace. Inverting the Cholesky factor costs
// another third of a cube and holds a second matrix the size of C, so a
// factorization that will only ever be solved against keeps the factor itself
// and substitutes through it instead. Set needTrace to false in that case.
class CoefMatFactor {
public:
  void compute(const arma::mat& C, bool needTrace = true){
    hasTrace = needTrace;
    // Factored in place, so the factor is never held twice
    useChol = arma::chol(Lfac, C, "lower");
    if(useChol && needTrace){
      useChol = arma::inv(Linv, arma::trimatl(Lfac));
      Lfac.reset();
    }
    if(useChol){
      Full.reset();
    }else{
      Linv.reset();
      Lfac.reset();
      hasTrace = true; // The fallback inverse supplies a trace either way
      Full = arma::inv(C);
    }
  }

  arma::mat solve(const arma::mat& b) const {
    if(useChol){
      if(hasTrace){
        return Linv.t()*(Linv*b);
      }
      // Forward substitution through L, then back substitution through L'.
      // Transposing the triangular view rather than the matrix keeps this
      // from building a transposed copy of the factor.
      return arma::solve(arma::trimatl(Lfac).t(),
                         arma::solve(arma::trimatl(Lfac), b));
    }
    return Full*b;
  }

  // Trace of the diagonal block of C^-1 spanning columns [first,last].
  // Only valid when the factorization was computed with needTrace.
  double blockTrace(arma::uword first, arma::uword last) const {
    if(useChol){
      return arma::accu(arma::square(Linv.cols(first,last)));
    }
    return arma::sum(Full(arma::span(first,last),
                          arma::span(first,last)).diag());
  }

private:
  arma::mat Linv; // inverse of the lower Cholesky factor, for a trace
  arma::mat Lfac; // the lower Cholesky factor itself, when no trace is wanted
  arma::mat Full; // full inverse, used only for the fallback
  bool useChol = false;
  bool hasTrace = true;
};

// Replacement for Armadillo's eig_sym
// Fixes an error with decomposition of large matrices
// If calcVec = false, eigvec is not used
// It would be better to template this function in the future
int eigen2(arma::vec& eigval, arma::mat& eigvec, arma::mat X,
           bool calcVec){
  char JOBZ;
  if(calcVec){
    JOBZ = 'V';
  }else{
    JOBZ = 'N';
  }
  char RANGE = 'A';
  char UPLO = 'L';
  long long int N = X.n_rows;
  // A = X
  long long int LDA = N;
  double VL = 0.0;
  double VU = 0.0;
  long long int IL = 0;
  long long int IU = 0;
  double ABSTOL = 0.0;
  long long int M = N;
  // W=eigval
  // Z=eigvec
  long long int LDZ = N;
  arma::Col<long long int> ISUPPZ(2*M);
  // WORK length to be determined
  double tmpWORK;
  long long int LWORK = -1; // To be calculated
  // IWORK length to be determined
  long long int tmpIWORK = 0;
  long long int LIWORK = -1; // To be calculated
  long long int INFO = 0;
  // Calculate LWORK and LIWORK
  F77_CALL(dsyevr)(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
          &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&tmpWORK,&LWORK,&tmpIWORK,&LIWORK,&INFO);
  LWORK = (long long int) tmpWORK;
  LIWORK = tmpIWORK;
  // Allocate WORK and IWORK
  arma::vec WORK(LWORK);
  arma::Col<long long int> IWORK(LIWORK);
  // Perform decomposition
  F77_CALL(dsyevr)(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
          &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&*WORK.begin(),&LWORK,&*IWORK.begin(),&LIWORK,&INFO);
  return INFO; // Return error code
}

// Objective function for REML using the FaST-LMM algorithm of Lippert et al.
// (2011), doi:10.1038/nmeth.1681. Written from the published description of
// the method; no code is taken from the FaST-LMM software.
//
// Unlike objREML and objREML2, the fixed effects are not projected out of
// the problem in advance. They are re-estimated by generalized least squares
// at every delta, which is why the restricted likelihood carries the
// log|X'*H^-1*X| term explicitly. The two formulations differ by a constant
// and are minimized at the same delta.
//
// param is log(delta). The parameter spans twenty decades and the likelihood
// is flat over most of them on a linear scale, so the search is made on the
// log scale.
Rcpp::List objREMLFastLMM(double param, Rcpp::List args){
  double df = args["df"];
  arma::vec yt = args["yt"];
  arma::mat Xt = args["Xt"];
  arma::vec s2 = args["s2"];
  double yy0 = args["yy0"];
  arma::vec Xy0 = args["Xy0"];
  arma::mat XX0 = args["XX0"];
  double nNull = args["nNull"];
  FastLMMFit fit = fastLMMFit(exp(param), yt, Xt, s2, yy0, Xy0, XX0, nNull);
  double value;
  if(fit.ok){
    value = df*log(fit.r) + fit.logDetH + fit.logDetA;
  }else{
    // X is rank deficient at this delta, so the point is excluded
    value = 1.0e100;
  }
  return Rcpp::List::create(Rcpp::Named("objective") = value,
                            Rcpp::Named("output") = 0);
}

// Finds the delta that minimizes the FaST-LMM REML objective.
//
// A coarse grid over log(delta) brackets the optimum before the bracket is
// refined, as recommended for FaST-LMM. Each grid point is linear in the rank
// of M*M', so the grid costs far less than the decomposition that precedes
// it, and it guards against the search settling on a local optimum.
double fastLMMDelta(Rcpp::List args){
  // delta = Ve/Vu carries the scale of the kernel, which is arbitrary:
  // multiplying the kernel by c multiplies the optimal delta by c and leaves
  // the fit unchanged. The search therefore runs from ten decades below the
  // smallest positive eigenvalue to ten above the largest, rather than around
  // one. A single average is not enough: a large eigenvalue whose direction
  // is absorbed by the fixed effects would set the scale while the variance
  // is estimated from the small ones, and the optimum could fall outside the
  // range. Eigenvalues within rounding of zero are not positive for this
  // purpose, since they would only widen the range.
  arma::vec s2 = args["s2"];
  double sMax = (s2.n_elem>0) ? s2.max() : 0.0;
  double sMin = sMax;
  double sTol = double(s2.n_elem)*arma::datum::eps*sMax;
  for(arma::uword i=0; i<s2.n_elem; ++i){
    if((s2(i)>sTol) && (s2(i)<sMin)) sMin = s2(i);
  }
  if(!(sMax>0.0) || !std::isfinite(sMax) || !(sMin>0.0)){
    sMin = 1.0;
    sMax = 1.0;
  }
  double logLower = log(1.0e-10)+log(sMin);
  double logUpper = log(1.0e10)+log(sMax);
  // Five grid points per decade, which is the density of the original
  // twenty-decade grid of 100 points
  int nGrid = std::max(100, int(std::ceil(5.0*(logUpper-logLower)/log(10.0)))+1);
  double step = (logUpper-logLower)/double(nGrid-1);
  int best = 0;
  double bestVal = 0.0;
  for(int i=0; i<nGrid; ++i){
    Rcpp::List gridOut = objREMLFastLMM(logLower+double(i)*step, args);
    double value = gridOut["objective"];
    if( (i==0) || (value<bestVal) ){
      bestVal = value;
      best = i;
    }
  }
  int lowIdx = (best>0) ? (best-1) : 0;
  int highIdx = (best<(nGrid-1)) ? (best+1) : (nGrid-1);
  Rcpp::List optRes = optimize(*objREMLFastLMM, args,
                               logLower+double(lowIdx)*step,
                               logLower+double(highIdx)*step,
                               1000, false, true, true);
  double logDelta = optRes["parameter"];
  return exp(logDelta);
}

// Produces a sum to zero design matrix with an intercept
arma::mat makeX(arma::uvec& x){
  arma::uword nTrain = x.n_elem;
  arma::uword nLevels = x.max();
  arma::mat X(nTrain,nLevels);
  if(nLevels==1){
    X.ones();
  }else{
    X.zeros();
    X.col(0).ones();
    for(arma::uword i=0; i<nTrain; ++i){
      if(x(i)==nLevels){
        X(i,arma::span(1,nLevels-1)).fill(-1.0);
      }else{
        X(i,x(i))=1;
      }
    }
  }
  return X;
}

// Produces a genotype design matrix
// z is an indicator vector matching y to individuals in G
// nGeno is the number of genotypes in the G matrix
arma::mat makeZ(arma::uvec& z, arma::uword nGeno){
  arma::uword nTrain = z.n_elem;
  arma::mat Z(nTrain,nGeno,arma::fill::zeros);
  for(arma::uword i=0; i<nTrain; ++i){
    Z(i,z(i)) = 1;
  }
  return Z;
}


//' @title Solve RR-BLUP
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e}, where
//' \eqn{u \sim N(0,I\sigma^2_u)} and \eqn{e \sim N(0,I\sigma^2_e)}. The
//' variance components are estimated by REML
//' \insertCite{patterson1971Recovery}{AlphaSimR} using the factored spectral
//' approach of FaST-LMM \insertCite{lippert2011FaST}{AlphaSimR}.
//'
//' @details
//' REML is reduced to a one dimensional search over delta, the ratio of the
//' residual variance to the marker variance. FaST-LMM works with M*M' itself
//' and re-estimates the fixed effects at every delta. One decomposition
//' therefore supplies the likelihood, the generalized least squares solution
//' for the fixed effects, and the BLUPs, and that decomposition is taken in
//' whichever of the two spaces is smaller. When there are fewer markers than
//' records the eigenvectors of M'*M serve in place of those of M*M', the rank
//' deficient directions are summarized analytically, and no matrix larger
//' than M is ever formed.
//'
//' This is an independent implementation of the published method. It shares
//' no code with the FaST-LMM software distributed by Microsoft.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns. Columns that are linearly
//' dependent on the others are dropped and get a fixed effect of zero. An
//' intercept, a column whose values are all the same, is always kept.
//' @param M a matrix with n rows and m columns
//'
//' @return A list with the marker variance \code{Vu}, the residual variance
//' \code{Ve}, the fixed effects \code{beta} and the marker effects \code{u}.
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M){
  arma::uword n = y.n_rows;
  arma::uword m = M.n_cols;
  if( (X.n_rows!=n) || (M.n_rows!=n) ){
    Rcpp::stop("solveRRBLUP: y, X and M must have the same number of rows");
  }
  if(!y.is_finite() || !M.is_finite()){
    Rcpp::stop("solveRRBLUP: y and M must not contain missing or infinite "
               "values");
  }
  // Columns of X that depend on the others are dropped, and the fit uses Xr,
  // an orthonormal basis for the span of the rest, so a badly scaled X costs
  // no precision; see fixedDesign
  const FixedDesign fd = fixedDesign(X);
  if(fd.keep.n_elem==0){
    Rcpp::stop("solveRRBLUP: X must have at least one nonzero column");
  }
  const arma::mat& Xr = fd.Q;
  arma::uword q = Xr.n_cols;
  if(q>=n){
    Rcpp::stop("solveRRBLUP: the rank of X must be smaller than the number "
               "of records");
  }
  double df = double(n)-double(q);

  // A least squares fit is taken off y and added back to the fixed effects
  // at the end, which keeps a large mean from swamping the sums of squares
  arma::mat b0;
  arma::vec yv = removeLeastSquares(Xr, y, b0).col(0);

  arma::vec s2;  // eigenvalues of M*M' that are kept
  arma::vec yt;  // U'*y
  arma::mat Xt;  // U'*X
  arma::mat MtU; // M'*U, low rank branch only
  arma::mat U;   // full rank branch only
  // Quantities outside the span of the kept eigenvectors
  double yy0 = 0.0;
  arma::vec Xy0(q, arma::fill::zeros);
  arma::mat XX0(q, q, arma::fill::zeros);
  double nNull = 0.0;

  if(m<n){
    // M'*M and M*M' share their nonzero eigenvalues, and the eigenvectors of
    // the larger follow from those of the smaller as U = M*V*diag(1/s). The
    // rotations below use that identity, so U itself is never formed.
    arma::vec Mty = M.t()*yv;
    arma::mat MtX = M.t()*Xr;
    arma::vec eigval(m);
    arma::mat eigvec(m,m);
    eigen2(eigval, eigvec, M.t()*M);
    double cutoff = rankCutoff(eigval.max(), n);
    arma::uvec keep = find(eigval>cutoff);
    s2 = eigval(keep);
    arma::mat Vk = eigvec.cols(keep);
    arma::vec sInv = 1.0/sqrt(s2);
    yt = sInv%(Vk.t()*Mty);
    Xt = Vk.t()*MtX;
    Xt.each_col() %= sInv;
    MtU = Vk;
    MtU.each_row() %= sqrt(s2).t();
    // The remaining n-k directions have a zero eigenvalue. Their
    // contribution is the same for every delta apart from a factor of
    // 1/delta, so it is summarized here rather than being stored.
    nNull = double(n)-double(s2.n_elem);
    // The parts of y and X in the null space are formed before their cross
    // products. Subtracting cross products instead, as in X'*X-Xt'*Xt, can
    // leave XX0 slightly negative when X lies in the span of the kernel,
    // and a tiny delta then turns that into a negative precision for the
    // fixed effects. The projection on the span is applied through the
    // rotated values already in hand, U*U'*v = M*Vk*diag(1/sqrt(s2))*(U'*v),
    // where U'*v is yt or Xt.
    arma::vec y0 = yv-M*(Vk*(sInv%yt));
    arma::mat X0 = Xr-M*(Vk*(Xt.each_col()%sInv));
    yy0 = arma::dot(y0, y0);
    Xy0 = X0.t()*y0;
    XX0 = X0.t()*X0;
  }else{
    // M*M' is the smaller of the two, and its decomposition spans all of
    // R^n, so there is no null space to account for
    arma::mat K = M*M.t();
    arma::vec eigval(n);
    U.set_size(n,n);
    eigen2(eigval, U, K);
    K.reset();
    eigval.elem(find(eigval<=rankCutoff(eigval.max(), n))).zeros();
    s2 = eigval;
    yt = U.t()*yv;
    Xt = U.t()*Xr;
  }

  Rcpp::List args = Rcpp::List::create(Rcpp::Named("df")=df,
                                       Rcpp::Named("yt")=yt,
                                       Rcpp::Named("Xt")=Xt,
                                       Rcpp::Named("s2")=s2,
                                       Rcpp::Named("yy0")=yy0,
                                       Rcpp::Named("Xy0")=Xy0,
                                       Rcpp::Named("XX0")=XX0,
                                       Rcpp::Named("nNull")=nNull);

  // Estimate variances. A coarse grid over log(delta) brackets the optimum
  // before the bracket is refined, as recommended for FaST-LMM. Each grid
  // point is linear in the rank of M*M', so the grid costs far less than the
  // decomposition that precedes it, and it guards against the search
  // settling on a local optimum.
  double delta = fastLMMDelta(args);

  // Solve the mixed model equations at the estimated delta
  FastLMMFit fit = fastLMMFit(delta, yt, Xt, s2, yy0, Xy0, XX0, nNull);
  if(!fit.ok){
    Rcpp::stop("solveRRBLUP: mixed model equations could not be solved");
  }
  arma::mat beta = fixedEffectsFromBasis(fd, fit.beta+b0, X.n_cols);

  // The BLUPs are u = M'*H^-1*(y-X*beta), which the decomposition already
  // in hand applies without a second factorization
  arma::vec et = yt-Xt*fit.beta; // U'*(y-X*beta)
  arma::vec etd = et/(s2+delta);
  arma::mat u;
  if(m<n){
    // u is MtU*etd, the sqrt(s2) in MtU cancelling the 1/sqrt(s2) of the
    // rotation. Eigenvectors of M'*M below the cutoff are given no effect:
    // the likelihood was fitted with their eigenvalues set to zero, so the
    // BLUPs are those of that same model. Along the kept ones the null
    // space term M'*(I-U*U')*e is exactly zero, and forming it as
    // M'*e-MtU*et would only leave rounding that 1/delta amplifies, which
    // ruined the effects when delta was tiny, as for a response with almost
    // no residual variation. Along the discarded ones 1/delta would also
    // stand in for 1/(s+delta), which is far too large when delta is tiny.
    u = MtU*etd;
  }else{
    // An eigenvalue set to zero by the rank rule was fitted as having no
    // genetic variance, so its direction gets no effect. Left in, it would
    // carry a gain of 1/delta where the kernel itself implies 1/(s+delta).
    etd.elem(find(s2==0.0)).zeros();
    u = M.t()*(U*etd);
  }

  double Vu = fit.r/df;
  double Ve = delta*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u);
}

//' @title Solve RR-BLUP with EM
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e} using
//' the Expectation-Maximization algorithm.
//'
//' @param Y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param Vu initial guess for variance of marker effects
//' @param Ve initial guess for error variance
//' @param tol tolerance for declaring convergence
//' @param maxIter maximum iteration for attempting convergence
//' @param useEM should EM algorithm be used. If false, no estimation of
//' variance components is performed. The initial values are treated as true.
//'
//' @details
//' The model is solved in one of two ways. Henderson's mixed model equations
//' work with a square matrix of the fixed effects plus the markers. Written
//' on the records instead, the model works with a square matrix of the
//' records. The two give the same answer, and the cost of each is the cost of
//' the other with the number of records and the number of markers exchanged,
//' so the smaller of those two numbers names the cheaper method.
//'
//' Estimating variance components requires the trace of the marker block of
//' the inverse coefficient matrix, which only Henderson's equations supply, so
//' \code{useEM = TRUE} always uses the marker side however many records there
//' are. With \code{useEM = FALSE} the choice is free and the cheaper method is
//' taken.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP_EM(arma::mat& Y, arma::mat& X,
                          arma::mat& M, double Vu, double Ve, 
                          double tol, int maxIter,
                          bool useEM){
  double lambda = Ve/Vu;
  double delta=0,VeN=0,VuN=0;
  int iter=0;
  arma::uword n=Y.n_rows,m=M.n_cols,q=X.n_cols;
  // The EM algorithm needs the trace of the marker block of the inverse
  // coefficient matrix, which is a property of Henderson's equations, so
  // estimating variance components always uses the marker side. With the
  // variance components held fixed either method gives the same answer and
  // the cheaper one is taken. See useRecordSide for how that is decided.
  if(!useEM && useRecordSide(n,m)){
    // V = M*M'*Vu + I*Ve, factored once and solved against the fixed effects
    // and the response together. Nothing of size m by n is formed: M'
    // multiplies a single vector of records rather than the whole of V^-1.
    arma::mat V = M*M.t();
    V *= Vu;
    V.diag() += Ve;
    arma::mat rhs = arma::join_rows(X,Y);
    arma::mat Vinvrhs;
    if(!arma::solve(Vinvrhs, V, rhs, arma::solve_opts::likely_sympd)){
      Rcpp::stop("solveRRBLUP_EM: mixed model equations could not be solved");
    }
    V.reset();
    rhs.reset();
    const arma::mat VinvX = Vinvrhs.head_cols(q);
    const arma::mat VinvY = Vinvrhs.tail_cols(1);
    arma::mat beta = solve(X.t()*VinvX, X.t()*VinvY);
    arma::mat u = M.t()*(VinvY-VinvX*beta)*Vu;
    return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                              Rcpp::Named("Ve")=Ve,
                              Rcpp::Named("beta")=beta,
                              Rcpp::Named("u")=u,
                              Rcpp::Named("iter")=iter);
  }
  arma::mat RHS(q+m,q+m),LHS(q+m,1),Rvec(q+m,1);
  const arma::mat XtM = X.t()*M;
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m-1)) = XtM;
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) = XtM.t();
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M.t()*M;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M.t()*Y;
  CoefMatFactor RHSfac;
  // Only the EM iterations ask for a trace, so a single solve keeps the
  // Cholesky factor rather than paying to invert it
  RHSfac.compute(RHS, useEM);
  LHS = RHSfac.solve(Rvec);
  if(useEM){
    const double YtY = as_scalar(Y.t()*Y);
    VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
    VuN = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*RHSfac.blockTrace(q, q+m-1)
    )/m;
    delta = VeN/VuN-lambda;
    while(fabs(delta)>tol){
      Ve = VeN;
      Vu = VuN;
      RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += delta;
      lambda += delta;
      RHSfac.compute(RHS);
      LHS = RHSfac.solve(Rvec);
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
      VuN = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*RHSfac.blockTrace(q, q+m-1)
      )/m;
      delta = VeN/VuN-lambda;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=LHS.rows(arma::span(0,q-1)),
                            Rcpp::Named("u")=LHS.rows(arma::span(q,q+m-1)),
                            Rcpp::Named("iter")=iter);
}

//' @title Solve RR-BLUP with EM and 2 random effects
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+M_1u_1+M_2u_2+e} using
//' the Expectation-Maximization algorithm.
//'
//' @param Y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M1 a matrix with n rows and m1 columns
//' @param M2 a matrix with n rows and m2 columns
//' @param Vu1 initial guess for variance of the first marker effects
//' @param Vu2 initial guess for variance of the second marker effects
//' @param Ve initial guess for error variance
//' @param tol tolerance for declaring convergence
//' @param maxIter maximum iteration for attempting convergence
//' @param useEM should EM algorithm be used. If false, no estimation of
//' variance components is performed. The initial values are treated as true.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP_EM2(const arma::mat& Y, const arma::mat& X,
                           const arma::mat& M1, const arma::mat& M2, 
                           double Vu1, double Vu2, double Ve, 
                           double tol, int maxIter, bool useEM){
  double lambda1 = Ve/Vu1;
  double lambda2 = Ve/Vu2;
  double delta1=0,delta2=0,VeN=0,Vu1N=0,Vu2N=0;
  int iter=0;
  arma::uword n=Y.n_rows,m=M1.n_cols,q=X.n_cols;
  if(!useEM & (n<(2*m))){
    arma::mat Vinv = inv_sympd(M1*M1.t()*Vu1+M2*M2.t()*Vu2+arma::eye(n,n)*Ve);
    arma::mat beta = solve(X.t()*Vinv*X, X.t()*Vinv*Y);
    arma::mat u(m,2);
    u.col(0) = M1.t()*Vinv*(Y-X*beta)*Vu1;
    u.col(1) = M2.t()*Vinv*(Y-X*beta)*Vu2;
    arma::vec Vu(2);
    Vu(0) = Vu1;
    Vu(1) = Vu2;
    return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                              Rcpp::Named("Ve")=Ve,
                              Rcpp::Named("beta")=beta,
                              Rcpp::Named("u")=u,
                              Rcpp::Named("iter")=iter);
  }
  arma::mat RHS(q+2*m,q+2*m),LHS(q+2*m,1),Rvec(q+2*m,1);
  // Top row
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m-1)) = X.t()*M1;
  RHS(arma::span(0,q-1),arma::span(q+m,q+2*m-1)) = X.t()*M2;
  // Second row
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) =
    RHS(arma::span(0,q-1),arma::span(q,q+m-1)).t();
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M1.t()*M1;
  RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)) = M1.t()*M2;
  // Third row
  RHS(arma::span(q+m,q+2*m-1),arma::span(0,q-1)) =
    RHS(arma::span(0,q-1),arma::span(q+m,q+2*m-1)).t();
  RHS(arma::span(q+m,q+2*m-1),arma::span(q,q+m-1)) =
    RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)).t();
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)) = M2.t()*M2;
  // Add to diagonal
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += lambda2;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M1.t()*Y;
  Rvec(arma::span(q+m,q+2*m-1),0) = M2.t()*Y;
  CoefMatFactor RHSfac;
  // Only the EM iterations ask for a trace, so a single solve keeps the
  // Cholesky factor rather than paying to invert it
  RHSfac.compute(RHS, useEM);
  LHS = RHSfac.solve(Rvec);
  if(useEM){
    const double YtY = as_scalar(Y.t()*Y);
    VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
    Vu1N = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*RHSfac.blockTrace(q, q+m-1)
    )/m;
    Vu2N = as_scalar(
      LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
        Ve*RHSfac.blockTrace(q+m, q+2*m-1)
    )/m;
    delta1 = VeN/Vu1N-lambda1;
    delta2 = VeN/Vu2N-lambda2;
    while((fabs(delta1)>tol) || (fabs(delta2)>tol)){
      Ve = VeN;
      Vu1 = Vu1N;
      Vu2 = Vu2N;
      RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += delta1;
      RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += delta2;
      lambda1 += delta1;
      lambda2 += delta2;
      RHSfac.compute(RHS);
      LHS = RHSfac.solve(Rvec);
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
      Vu1N = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*RHSfac.blockTrace(q, q+m-1)
      )/m;
      Vu2N = as_scalar(
        LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
          Ve*RHSfac.blockTrace(q+m, q+2*m-1)
      )/m;
      delta1 = VeN/Vu1N-lambda1;
      delta2 = VeN/Vu2N-lambda2;
    }
  }
  arma::vec Vu(2);
  Vu(0) = Vu1;
  Vu(1) = Vu2;
  arma::mat u(m,2);
  u.col(0) = LHS.rows(arma::span(q,q+m-1));
  u.col(1) = LHS.rows(arma::span(q+m,q+2*m-1));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=LHS.rows(arma::span(0,q-1)),
                            Rcpp::Named("u")=u,
                            Rcpp::Named("iter")=iter);
}

//' @title Solve RR-BLUP with EM and 3 random effects
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+M_1u_1+M_2u_2+M_3u_3+e} using
//' the Expectation-Maximization algorithm.
//'
//' @param Y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M1 a matrix with n rows and m1 columns
//' @param M2 a matrix with n rows and m2 columns
//' @param M3 a matrix with n rows and m3 columns
//' @param Vu1 initial guess for variance of the first marker effects
//' @param Vu2 initial guess for variance of the second marker effects
//' @param Vu3 initial guess for variance of the third marker effects
//' @param Ve initial guess for error variance
//' @param tol tolerance for declaring convergence
//' @param maxIter maximum iteration for attempting convergence
//' @param useEM should EM algorithm be used. If false, no estimation of
//' variance components is performed. The initial values are treated as true.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP_EM3(const arma::mat& Y, const arma::mat& X,
                           const arma::mat& M1, const arma::mat& M2, 
                           const arma::mat& M3, double Vu1, double Vu2, 
                           double Vu3, double Ve, double tol, 
                           int maxIter, bool useEM){
  double lambda1 = Ve/Vu1;
  double lambda2 = Ve/Vu2;
  double lambda3 = Ve/Vu3;
  double delta1=0,delta2=0,delta3=0,VeN=0,Vu1N=0,Vu2N=0,Vu3N=0;
  int iter=0;
  arma::uword n=Y.n_rows,m=M1.n_cols,q=X.n_cols;
  if(!useEM & (n<(3*m))){
    arma::mat Vinv = inv_sympd(M1*M1.t()*Vu1+M2*M2.t()*Vu2+M3*M3.t()*Vu3+arma::eye(n,n)*Ve);
    arma::mat beta = solve(X.t()*Vinv*X, X.t()*Vinv*Y);
    arma::mat u(m,3);
    u.col(0) = M1.t()*Vinv*(Y-X*beta)*Vu1;
    u.col(1) = M2.t()*Vinv*(Y-X*beta)*Vu2;
    u.col(2) = M3.t()*Vinv*(Y-X*beta)*Vu3;
    arma::vec Vu(3);
    Vu(0) = Vu1;
    Vu(1) = Vu2;
    Vu(2) = Vu3;
    return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                              Rcpp::Named("Ve")=Ve,
                              Rcpp::Named("beta")=beta,
                              Rcpp::Named("u")=u,
                              Rcpp::Named("iter")=iter);
  }
  arma::mat RHS(q+3*m,q+3*m),LHS(q+3*m,1),Rvec(q+3*m,1);
  // Top row
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m-1)) = X.t()*M1;
  RHS(arma::span(0,q-1),arma::span(q+m,q+2*m-1)) = X.t()*M2;
  RHS(arma::span(0,q-1),arma::span(q+2*m,q+3*m-1)) = X.t()*M3;
  // Second row
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) =
    RHS(arma::span(0,q-1),arma::span(q,q+m-1)).t();
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M1.t()*M1;
  RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)) = M1.t()*M2;
  RHS(arma::span(q,q+m-1),arma::span(q+2*m,q+3*m-1)) = M1.t()*M3;
  // Third row
  RHS(arma::span(q+m,q+2*m-1),arma::span(0,q-1)) =
    RHS(arma::span(0,q-1),arma::span(q+m,q+2*m-1)).t();
  RHS(arma::span(q+m,q+2*m-1),arma::span(q,q+m-1)) =
    RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)).t();
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)) = M2.t()*M2;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+2*m,q+3*m-1)) = M2.t()*M3;
  // Fourth row
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(0,q-1)) =
    RHS(arma::span(0,q-1),arma::span(q+2*m,q+3*m-1)).t();
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q,q+m-1)) =
    RHS(arma::span(q,q+m-1),arma::span(q+2*m,q+3*m-1)).t();
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+m,q+2*m-1)) =
    RHS(arma::span(q+m,q+2*m-1),arma::span(q+2*m,q+3*m-1)).t();
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)) = M3.t()*M3;
  // Add to diagonal
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += lambda2;
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)).diag() += lambda3;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M1.t()*Y;
  Rvec(arma::span(q+m,q+2*m-1),0) = M2.t()*Y;
  Rvec(arma::span(q+2*m,q+3*m-1),0) = M3.t()*Y;
  CoefMatFactor RHSfac;
  // Only the EM iterations ask for a trace, so a single solve keeps the
  // Cholesky factor rather than paying to invert it
  RHSfac.compute(RHS, useEM);
  LHS = RHSfac.solve(Rvec);
  if(useEM){
    const double YtY = as_scalar(Y.t()*Y);
    VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
    Vu1N = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*RHSfac.blockTrace(q, q+m-1)
    )/m;
    Vu2N = as_scalar(
      LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
        Ve*RHSfac.blockTrace(q+m, q+2*m-1)
    )/m;
    Vu3N = as_scalar(
      LHS(arma::span(q+2*m,q+3*m-1),0).t()*LHS(arma::span(q+2*m,q+3*m-1),0)+
        Ve*RHSfac.blockTrace(q+2*m, q+3*m-1)
    )/m;
    delta1 = VeN/Vu1N-lambda1;
    delta2 = VeN/Vu2N-lambda2;
    delta3 = VeN/Vu3N-lambda3;
    while((fabs(delta1)>tol) || (fabs(delta2)>tol) || (fabs(delta3)>tol)){
      Ve = VeN;
      Vu1 = Vu1N;
      Vu2 = Vu2N;
      Vu3 = Vu3N;
      RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += delta1;
      RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += delta2;
      RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)).diag() += delta3;
      lambda1 += delta1;
      lambda2 += delta2;
      lambda3 += delta3;
      RHSfac.compute(RHS);
      LHS = RHSfac.solve(Rvec);
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(YtY-LHS.t()*Rvec)/(n-q);
      Vu1N = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*RHSfac.blockTrace(q, q+m-1)
      )/m;
      Vu2N = as_scalar(
        LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
          Ve*RHSfac.blockTrace(q+m, q+2*m-1)
      )/m;
      Vu3N = as_scalar(
        LHS(arma::span(q+2*m,q+3*m-1),0).t()*LHS(arma::span(q+2*m,q+3*m-1),0)+
          Ve*RHSfac.blockTrace(q+2*m, q+3*m-1)
      )/m;
      delta1 = VeN/Vu1N-lambda1;
      delta2 = VeN/Vu2N-lambda2;
      delta3 = VeN/Vu3N-lambda3;
    }
  }
  arma::vec Vu(3);
  Vu(0) = Vu1;
  Vu(1) = Vu2;
  Vu(2) = Vu3;
  arma::mat u(m,3);
  u.col(0) = LHS.rows(arma::span(q,q+m-1));
  u.col(1) = LHS.rows(arma::span(q+m,q+2*m-1));
  u.col(2) = LHS.rows(arma::span(q+2*m,q+3*m-1));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=LHS.rows(arma::span(0,q-1)),
                            Rcpp::Named("u")=u,
                            Rcpp::Named("iter")=iter);
}

// The conjugate gradient iteration reads the marker genotypes twice per
// iteration, and at one byte per marker those reads are limited by memory
// bandwidth. The iteration therefore works from a copy packed to one bit per
// allele: marker j of the chip is bit j%8 of column j/8, and row i*P+h holds
// haplotype h of individual i, so an individual's haplotypes sit next to each
// other. For a diploid this is a quarter of the size of the dosage matrix.
// Bits past the last marker are zero.
//
// The copy is built a tile of individuals at a time. Markers that are close
// together in the genome then share the cache lines read for the tile, and
// each block writes only its own rows, so the result does not depend on the
// number of threads.
arma::Mat<unsigned char> packChipPlanes(const PackedGeno& g, const LociPos& pos,
                                        int nThreads){
  arma::uword ploidy = g.ploidy;
  arma::uword nInd = g.nInd;
  arma::uword nLoci = pos.base.size();
  arma::uword nCol = (nLoci+7)/8;
  arma::uword nRow = nInd*ploidy;
  arma::Mat<unsigned char> planes(nRow,nCol);
  const arma::uword kPackTile = 64;
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword t0=indStart; t0<indEnd; t0+=kPackTile){
      arma::uword t1 = std::min(t0+kPackTile, indEnd);
      for(arma::uword c=0; c<nCol; ++c){
        arma::uword jEnd = std::min(8*c+8, nLoci);
        unsigned char* col = planes.colptr(c);
        for(arma::uword i=t0; i<t1; ++i){
          for(arma::uword h=0; h<ploidy; ++h){
            unsigned int b = 0;
            for(arma::uword j=8*c; j<jEnd; ++j){
              const unsigned char* p = pos.base[j]+pos.indStride(j)*i+
                pos.hapStride(j)*h;
              b |= ((unsigned int) ((*p>>pos.bit(j)) & 1u))<<(j-8*c);
            }
            col[i*ploidy+h] = (unsigned char) b;
          }
        }
      }
    }
  }
  return planes;
}

// Applies the centered genotype matrix to a vector, t = Mc*p, from the
// packed copy of the genotypes.
//
// The additive coding of a dosage x at ploidy P is (2/P)*x-1, which is a sum
// over haplotypes of 2/P per allele carried, less one. Centering takes Mmean
// off as well, so Mc*p is the sum of (2/P)*p_j over the alleles an individual
// carries, less the constant sum of (1+Mmean_j)*p_j, which is passed in as
// offset. Each half of a column covers four markers, and a table over the
// sixteen bit patterns of those four markers turns one haplotype's half byte
// into its contribution with one lookup, as gvByTable does for genetic
// values. The tables cost one addition per entry and are rebuilt for every
// product.
//
// The work is split over blocks of records. A block owns its own part of t
// and every record sums the columns in order, so nothing is accumulated
// across threads and the result does not depend on how many threads are
// used.
void multMc(arma::vec& t, const arma::Mat<unsigned char>& planes,
            arma::uword ploidy, const arma::vec& offset,
            const arma::vec& p, arma::mat& tabs, int nThreads){
  arma::uword nRow = planes.n_rows;
  arma::uword nInd = nRow/ploidy;
  arma::uword nCol = planes.n_cols;
  arma::uword nLoci = p.n_elem;
  double alleleScale = 2.0/double(ploidy);
  const double* pPtr = p.memptr();
  double shift = 0.0;
  for(arma::uword j=0; j<nLoci; ++j){
    shift += offset(j)*pPtr[j];
  }
  // Each table entry is the entry with its lowest set bit cleared plus that
  // bit's value. The low half byte's table is rows 0 to 15 of a column and
  // the high half byte's is rows 16 to 31.
  for(arma::uword c=0; c<nCol; ++c){
    double v[8];
    for(arma::uword b=0; b<8; ++b){
      arma::uword j = 8*c+b;
      v[b] = (j<nLoci) ? alleleScale*pPtr[j] : 0.0;
    }
    double* tab = tabs.colptr(c);
    for(arma::uword half=0; half<2; ++half){
      double* tb = tab+16*half;
      const double* vb = v+4*half;
      tb[0] = 0.0;
      for(unsigned int q=1; q<16; ++q){
        unsigned int low = q & (0u-q);
        unsigned int bit = (low==1u) ? 0 : ((low==2u) ? 1 : ((low==4u) ? 2 : 3));
        tb[q] = tb[q & (q-1u)]+vb[bit];
      }
    }
  }
  double* tPtr = t.memptr();
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword i=indStart; i<indEnd; ++i){
      tPtr[i] = -shift;
    }
    for(arma::uword c=0; c<nCol; ++c){
      const unsigned char* col = planes.colptr(c);
      const double* tab = tabs.colptr(c);
      if(ploidy==2){
        // The common case is written out so that it has no haplotype loop
        for(arma::uword i=indStart; i<indEnd; ++i){
          unsigned int b0 = col[2*i], b1 = col[2*i+1];
          tPtr[i] += (tab[b0 & 15u]+tab[16+(b0>>4)])+
            (tab[b1 & 15u]+tab[16+(b1>>4)]);
        }
      }else{
        for(arma::uword i=indStart; i<indEnd; ++i){
          double s = 0.0;
          for(arma::uword h=0; h<ploidy; ++h){
            unsigned int b = col[i*ploidy+h];
            s += tab[b & 15u]+tab[16+(b>>4)];
          }
          tPtr[i] += s;
        }
      }
    }
  }
}

// Applies the transpose of the centered genotype matrix, w = Mc'*t, from the
// packed copy of the genotypes.
//
// For each half of a column, the values of t are added into sixteen bins
// chosen by each haplotype's bit pattern, and a marker's sum is then the sum
// of the eight bins whose pattern has its bit set. Records alternate between
// two sets of bins, and in a diploid so do the two haplotypes, so that
// consecutive additions rarely wait on the same bin. The sets are merged in a
// fixed order. The coding and centering are applied as in multMc.
//
// Split over blocks of columns. A block works out its own entries of w from
// its own columns, so again nothing is accumulated across threads.
void multMcT(arma::vec& w, const arma::Mat<unsigned char>& planes,
             arma::uword ploidy, const arma::vec& offset,
             const arma::vec& t, int nThreads){
  arma::uword nRow = planes.n_rows;
  arma::uword nInd = nRow/ploidy;
  arma::uword nCol = planes.n_cols;
  arma::uword nLoci = w.n_elem;
  double alleleScale = 2.0/double(ploidy);
  const double* tPtr = t.memptr();
  double tSum = 0.0;
  for(arma::uword i=0; i<nInd; ++i){
    tSum += tPtr[i];
  }
  double* wPtr = w.memptr();
  arma::uword nBlocks = countBlocks(nCol);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword colStart = blockStart(nCol, nBlocks, block);
    arma::uword colEnd = blockStart(nCol, nBlocks, block+1);
    // Four sets of 32 bins: the low half byte's bins, then the high's
    double bins[4][32];
    for(arma::uword c=colStart; c<colEnd; ++c){
      const unsigned char* col = planes.colptr(c);
      for(arma::uword s=0; s<4; ++s){
        for(arma::uword q=0; q<32; ++q){
          bins[s][q] = 0.0;
        }
      }
      if(ploidy==2){
        arma::uword i = 0;
        for(; i+1<nInd; i+=2){
          double t0 = tPtr[i], t1 = tPtr[i+1];
          unsigned int a0 = col[2*i], a1 = col[2*i+1];
          unsigned int b0 = col[2*i+2], b1 = col[2*i+3];
          bins[0][a0 & 15u] += t0;
          bins[0][16+(a0>>4)] += t0;
          bins[1][a1 & 15u] += t0;
          bins[1][16+(a1>>4)] += t0;
          bins[2][b0 & 15u] += t1;
          bins[2][16+(b0>>4)] += t1;
          bins[3][b1 & 15u] += t1;
          bins[3][16+(b1>>4)] += t1;
        }
        if(i<nInd){
          double t0 = tPtr[i];
          unsigned int a0 = col[2*i], a1 = col[2*i+1];
          bins[0][a0 & 15u] += t0;
          bins[0][16+(a0>>4)] += t0;
          bins[1][a1 & 15u] += t0;
          bins[1][16+(a1>>4)] += t0;
        }
      }else{
        for(arma::uword i=0; i<nInd; ++i){
          double ti = tPtr[i];
          double* bi = bins[2*(i & 1u)];
          for(arma::uword h=0; h<ploidy; ++h){
            unsigned int b = col[i*ploidy+h];
            bi[b & 15u] += ti;
            bi[16+(b>>4)] += ti;
          }
        }
      }
      double merged[32];
      for(arma::uword q=0; q<32; ++q){
        merged[q] = (bins[0][q]+bins[1][q])+(bins[2][q]+bins[3][q]);
      }
      for(arma::uword b=0; b<8; ++b){
        arma::uword j = 8*c+b;
        if(j>=nLoci){
          break;
        }
        const double* mb = merged+16*(b/4);
        unsigned int bit = 1u<<(b%4);
        double sum = 0.0;
        for(unsigned int q=1; q<16; ++q){
          if(q & bit){
            sum += mb[q];
          }
        }
        wPtr[j] = alleleScale*sum-offset(j)*tSum;
      }
    }
  }
}

// Builds the genotype cross product K = Msub*Msub' for a subset of records.
//
// The genotypes are read a panel of markers at a time and turned into dosages
// as they are read, so the only doubles held are one panel, whatever the
// number of markers. Armadillo sends the rank update to syrk.
arma::mat gramFromGeno(const arma::Mat<unsigned char>& M,
                       const arma::uvec& subset,
                       const double* xaPtr){
  arma::uword nSub = subset.n_elem;
  arma::uword nLoci = M.n_cols;
  const arma::uword panelWidth = 256;
  arma::mat K(nSub,nSub,arma::fill::zeros);
  for(arma::uword start=0; start<nLoci; start+=panelWidth){
    arma::uword stop = start+panelWidth;
    if(stop>nLoci){
      stop = nLoci;
    }
    arma::mat panel(nSub,stop-start);
    for(arma::uword j=start; j<stop; ++j){
      const unsigned char* Mcol = M.colptr(j);
      double* panelCol = panel.colptr(j-start);
      for(arma::uword i=0; i<nSub; ++i){
        panelCol[i] = xaPtr[Mcol[subset(i)]];
      }
    }
    K += panel*panel.t();
  }
  return K;
}

// Estimates the variance components by REML from a genotype cross product,
// using the same FaST-LMM machinery as solveRRBLUP. Only the variance
// components are worked out here; the marker effects are solved for
// afterward, on the full data set.
//
// The genotypes do not need centering first. The restricted likelihood depends
// on them only through S*M*M'*S for the projector S, and S*M is the same
// whether or not the columns of M have had their means taken off, so long as
// X holds an intercept, which makeX guarantees.
void remlFromGram(double& Vu, double& Ve, const arma::vec& yv,
                  const arma::mat& X, arma::mat& K){
  arma::uword n = yv.n_elem;
  arma::uword q = X.n_cols;
  double df = double(n)-double(q);
  arma::vec eigval(n);
  arma::mat U(n,n);
  eigen2(eigval, U, K);
  K.reset();
  eigval.elem(find(eigval<=rankCutoff(eigval.max(), n))).zeros();
  // A least squares fit is taken off y, which leaves REML unchanged and
  // keeps a large mean from swamping the sums of squares
  arma::mat b0;
  arma::vec yt = U.t()*removeLeastSquares(X, yv, b0).col(0);
  arma::mat Xt = U.t()*X;
  U.reset();
  // K spans all of R^n, so there is no null space to summarize
  double yy0 = 0.0;
  arma::vec Xy0(q,arma::fill::zeros);
  arma::mat XX0(q,q,arma::fill::zeros);
  double nNull = 0.0;
  Rcpp::List args = Rcpp::List::create(Rcpp::Named("df")=df,
                                       Rcpp::Named("yt")=yt,
                                       Rcpp::Named("Xt")=Xt,
                                       Rcpp::Named("s2")=eigval,
                                       Rcpp::Named("yy0")=yy0,
                                       Rcpp::Named("Xy0")=Xy0,
                                       Rcpp::Named("XX0")=XX0,
                                       Rcpp::Named("nNull")=nNull);
  double delta = fastLMMDelta(args);
  FastLMMFit fit = fastLMMFit(delta, yt, Xt, eigval, yy0, Xy0, XX0, nNull);
  if(!fit.ok){
    Rcpp::stop("fastRRBLUP: could not estimate the variance components");
  }
  Vu = fit.r/df;
  Ve = delta*Vu;
}

// Called by fastRRBLUP function
// Solves the mixed model equations of an RR-BLUP model by preconditioned
// conjugate gradient, iterating on the genotypes rather than on a stored
// coefficient matrix. See Stranden and Lidauer (1999),
// doi:10.3168/jds.S0022-0302(99)75535-9.
// x is an indicator vector for the fixed effect levels, as in callRRBLUP
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP(arma::vec y, arma::uvec x,
                          const Rcpp::S4& pop,
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                          double Vu, double Ve, arma::uword maxIter,
                          bool estVarComp, arma::uvec subset, int nThreads){
  // The marker dosages are decoded straight from the population's packed
  // genotypes. Taking the genotypes as an arma::field would first copy every
  // chromosome, including all the loci that are not markers.
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  DosageTiles tiles(g, pos);
  arma::uword ploidy = g.ploidy;
  arma::uword nInd = g.nInd;
  arma::uword nLoci = pos.base.size();
  arma::Mat<unsigned char> M(nInd,nLoci);
  {
    // A block of individuals is written as a tile whose stride is the
    // column length of M, so decode fills its rows of M in place
    arma::uword nBlocks = countBlocks(nInd);
    int decodeThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(decodeThreads)){
      decodeThreads = static_cast<int>(nBlocks);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(decodeThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword indStart = blockStart(nInd, nBlocks, block);
      arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
      tiles.decode(indStart, indEnd-indStart, nInd, M.memptr()+indStart);
    }
  }
  double dP = double(ploidy);
  // Sum to zero design matrix with an intercept, matching callRRBLUP.
  // A single fixed effect level gives a single column of ones.
  arma::mat Xfix = makeX(x);
  arma::uword q = Xfix.n_cols;
  
  // Genotypes are turned into dosages by a lookup as they are read, so the
  // value of each of the ploidy+1 genotypes is worked out once here
  arma::vec xa(ploidy+1);
  for(arma::uword i=0; i<xa.n_elem; ++i){
    xa(i) = (double(i)-dP/2.0)*(2.0/dP);
  }
  const double* xaPtr = xa.memptr();
  
  // The fixed effects are absorbed. Writing S for I-X*(X'X)^-1*X', the
  // equations reduce to (Mc'*S*Mc+lambda*I)*u = Mc'*S*y, and the fixed
  // effects follow from beta = (X'X)^-1*X'*(y-Mc*u) once u is known.
  arma::mat XtXi = inv_sympd(Xfix.t()*Xfix);
  arma::vec Sy = y-Xfix*(XtXi*(Xfix.t()*y));
  arma::rowvec Xsum = sum(Xfix,0); // X'*1
  
  // Variance components, when they have not been supplied.
  //
  // These are worked out by REML on a subset of the records, because the
  // decomposition REML needs costs the square of the number of records in
  // memory and their cube in time, neither of which this function can afford
  // on a full data set. Variance components are nuisance parameters here and
  // are estimated far more precisely than a breeding program needs, so a few
  // thousand records give values good enough to shrink with. The marker
  // effects are always solved for on every record.
  //
  // Whichever of the two cross products is smaller is the one decomposed. If
  // there are fewer markers than sampled records, solveRRBLUP is given the
  // sampled genotypes and takes the marker by marker route itself. Otherwise
  // the record by record cross product is built a panel at a time.
  if(estVarComp){
    subset -= 1; // R to C++
    arma::uword nSub = subset.n_elem;
    arma::vec ySub = y.elem(subset);
    arma::mat XSub = Xfix.rows(subset);
    if(nLoci<nSub){
      arma::mat MSub(nSub,nLoci);
      for(arma::uword j=0; j<nLoci; ++j){
        const unsigned char* Mcol = M.colptr(j);
        double* MSubCol = MSub.colptr(j);
        for(arma::uword i=0; i<nSub; ++i){
          MSubCol[i] = xaPtr[Mcol[subset(i)]];
        }
      }
      Rcpp::List varComp = solveRRBLUP(ySub, XSub, MSub);
      Vu = varComp["Vu"];
      Ve = varComp["Ve"];
    }else{
      arma::mat K = gramFromGeno(M, subset, xaPtr);
      remlFromGram(Vu, Ve, ySub, XSub, K);
    }
    if(!(Vu>0.0)){
      Rcpp::stop("fastRRBLUP: estimated marker variance is not positive");
    }
  }
  double lambda = Ve/Vu;
  
  // One pass over the genotypes gives the column means, the right hand side
  // and the diagonal of the coefficient matrix, which is what the iteration
  // is preconditioned with. Everything is accumulated into scalars, so no
  // column of dosages is ever formed and the pass splits over loci without
  // needing memory per thread.
  arma::rowvec Mmean(nLoci);
  arma::vec r(nLoci), precond(nLoci);
  double SySum = accu(Sy);
  const double* SyPtr = Sy.memptr();
  const double* XfixPtr = Xfix.memptr();
  {
    arma::uword nBlocks = countBlocks(nLoci);
    int setupThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(setupThreads)){
      setupThreads = static_cast<int>(nBlocks);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(setupThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      arma::vec XtM(q);
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        const unsigned char* Mcol = M.colptr(j);
        double s1 = 0.0, s2 = 0.0, sSy = 0.0;
        XtM.zeros();
        for(arma::uword i=0; i<nInd; ++i){
          double value = xaPtr[Mcol[i]];
          s1 += value;
          s2 += value*value;
          sSy += value*SyPtr[i];
          for(arma::uword k=0; k<q; ++k){
            XtM(k) += XfixPtr[i+k*nInd]*value;
          }
        }
        double mean = s1/double(nInd);
        Mmean(j) = mean;
        r(j) = sSy-mean*SySum;
        // Centering the column shifts X'*Mc and takes n*mean^2 off Mc'*Mc
        for(arma::uword k=0; k<q; ++k){
          XtM(k) -= mean*Xsum(k);
        }
        double diagonal = s2-s1*mean-
          as_scalar(XtM.t()*XtXi*XtM)+lambda;
        if(diagonal<lambda){
          // Rounding only, the true value cannot be below lambda
          diagonal = lambda;
        }
        precond(j) = diagonal;
      }
    }
  }
  
  // The dosages have done their job. The iteration works from a copy of the
  // markers packed to one bit per allele, built only once M has been freed,
  // so the most memory held at once is still M's.
  M.reset();
  arma::Mat<unsigned char> planes = packChipPlanes(g, pos, nThreads);
  arma::vec offset = 1.0+Mmean.t();
  arma::mat tabs(32,planes.n_cols);

  // Preconditioned conjugate gradient. Every iteration is one product with
  // the coefficient matrix, taken as Mc'*(S*(Mc*p))+lambda*p, so the matrix
  // itself is never formed. The vectors below are the only working memory
  // the solver needs.
  arma::vec u(nLoci,arma::fill::zeros);
  arma::vec z = r/precond;
  arma::vec pDir = z;
  arma::vec Ap(nLoci);
  arma::vec t(nInd);
  double rNorm0 = norm(r,2);
  double rz = dot(r,z);
  bool converged = (rNorm0<=0.0);
  for(arma::uword iter=0; (iter<maxIter)&&(!converged); ++iter){
    multMc(t, planes, ploidy, offset, pDir, tabs, nThreads);
    t -= Xfix*(XtXi*(Xfix.t()*t)); // S*t
    multMcT(Ap, planes, ploidy, offset, t, nThreads);
    Ap += lambda*pDir;
    double pAp = dot(pDir,Ap);
    if(pAp<=0.0){
      // The matrix is positive definite, so this is rounding once there is
      // nothing left to solve for
      converged = true;
      break;
    }
    double step = rz/pAp;
    u += step*pDir;
    r -= step*Ap;
    if(norm(r,2)<=1.0e-8*rNorm0){
      converged = true;
      break;
    }
    z = r/precond;
    double rzNew = dot(r,z);
    pDir = z+(rzNew/rz)*pDir;
    rz = rzNew;
  }
  if(!converged){
    Rf_warning("callFastRRBLUP: reached maxIter without converging");
  }
  
  // Fixed effects, recovered from the solution. They are found against the
  // centered genotypes, while the genetic values are predicted from the
  // uncentered ones, so the intercept takes off Mmean*u.
  multMc(t, planes, ploidy, offset, u, tabs, nThreads);
  arma::vec beta = XtXi*(Xfix.t()*(y-t));
  double shift = as_scalar(Mmean*u);
  
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-shift,
                            Rcpp::Named("mu")=beta(0)-shift,
                            Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve);
}

// The marker dosages of a chip, one byte per marker, decoded straight from the
// population's packed genotypes. A block of individuals is written as a tile
// whose stride is the column length of M, so decode fills its rows of M in
// place.
arma::Mat<unsigned char> decodeChipDosage(const PackedGeno& g,
                                          const LociPos& pos, int nThreads){
  DosageTiles tiles(g, pos);
  arma::uword nInd = g.nInd;
  arma::uword nLoci = pos.base.size();
  arma::Mat<unsigned char> M(nInd,nLoci);
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    tiles.decode(indStart, indEnd-indStart, nInd, M.memptr()+indStart);
  }
  return M;
}

// Applies the centered genotype matrix to several vectors at once,
// Tt = (Mc*Pt')', from the packed copy of the genotypes.
//
// This is multMc with one column of Pt per marker and one row per trait. The
// tables hold every trait's contribution for each bit pattern side by side,
// so the half byte of a haplotype is read and looked up once for all of the
// traits, and an individual's results are contiguous in Tt. The traits share
// the reads of the genotypes, but with the markers packed those reads are
// cheap, and the additions, which grow with the number of traits, take most
// of the time. A product for several traits therefore costs about as much as
// one product per trait.
//
// The work is split over blocks of records. A block owns its own columns of
// Tt and every record sums the columns of the genotypes in order, so the
// result does not depend on the number of threads.
void multMcBlock(arma::mat& Tt, const arma::Mat<unsigned char>& planes,
                 arma::uword ploidy, const arma::vec& offset,
                 const arma::mat& Pt, arma::mat& tabs, int nThreads){
  arma::uword nT = Pt.n_rows;
  arma::uword nLoci = Pt.n_cols;
  arma::uword nRow = planes.n_rows;
  arma::uword nInd = nRow/ploidy;
  arma::uword nCol = planes.n_cols;
  double alleleScale = 2.0/double(ploidy);
  arma::vec shift(nT, arma::fill::zeros);
  for(arma::uword j=0; j<nLoci; ++j){
    for(arma::uword k=0; k<nT; ++k){
      shift(k) += offset(j)*Pt(k,j);
    }
  }
  // As in multMc, each table entry is the entry with its lowest set bit
  // cleared plus that bit's value, here for every trait at once. Entry q of
  // the low half byte's table starts at q*nT and of the high one's at
  // (16+q)*nT.
  std::vector<double> v(8*nT);
  for(arma::uword c=0; c<nCol; ++c){
    for(arma::uword b=0; b<8; ++b){
      arma::uword j = 8*c+b;
      for(arma::uword k=0; k<nT; ++k){
        v[b*nT+k] = (j<nLoci) ? alleleScale*Pt(k,j) : 0.0;
      }
    }
    double* tab = tabs.colptr(c);
    for(arma::uword half=0; half<2; ++half){
      double* tb = tab+16*half*nT;
      const double* vb = v.data()+4*half*nT;
      for(arma::uword k=0; k<nT; ++k){
        tb[k] = 0.0;
      }
      for(unsigned int q=1; q<16; ++q){
        unsigned int low = q & (0u-q);
        unsigned int bit = (low==1u) ? 0 : ((low==2u) ? 1 : ((low==4u) ? 2 : 3));
        const double* prev = tb+(q & (q-1u))*nT;
        const double* add = vb+bit*nT;
        double* out = tb+q*nT;
        for(arma::uword k=0; k<nT; ++k){
          out[k] = prev[k]+add[k];
        }
      }
    }
  }
  arma::uword nBlocks = countBlocks(nInd);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword indStart = blockStart(nInd, nBlocks, block);
    arma::uword indEnd = blockStart(nInd, nBlocks, block+1);
    for(arma::uword i=indStart; i<indEnd; ++i){
      double* ti = Tt.colptr(i);
      for(arma::uword k=0; k<nT; ++k){
        ti[k] = -shift(k);
      }
    }
    for(arma::uword c=0; c<nCol; ++c){
      const unsigned char* col = planes.colptr(c);
      const double* tab = tabs.colptr(c);
      for(arma::uword i=indStart; i<indEnd; ++i){
        double* ti = Tt.colptr(i);
        for(arma::uword h=0; h<ploidy; ++h){
          unsigned int b = col[i*ploidy+h];
          const double* lo = tab+(b & 15u)*nT;
          const double* hi = tab+(16+(b>>4))*nT;
          for(arma::uword k=0; k<nT; ++k){
            ti[k] += lo[k]+hi[k];
          }
        }
      }
    }
  }
}

// Applies the transpose of the centered genotype matrix to several vectors at
// once, Wt = (Mc'*Tt')', from the packed copy of the genotypes.
//
// This is multMcT with every bin holding one value per trait. Records
// alternate between two pairs of bin sets, as there, and the sets are merged
// in a fixed order. Split over blocks of columns, each of which works out its
// own columns of Wt, so nothing is accumulated across threads.
void multMcTBlock(arma::mat& Wt, const arma::Mat<unsigned char>& planes,
                  arma::uword ploidy, const arma::vec& offset,
                  const arma::mat& Tt, int nThreads){
  arma::uword nT = Tt.n_rows;
  arma::uword nInd = Tt.n_cols;
  arma::uword nCol = planes.n_cols;
  arma::uword nLoci = Wt.n_cols;
  double alleleScale = 2.0/double(ploidy);
  arma::vec tSum(nT, arma::fill::zeros);
  for(arma::uword i=0; i<nInd; ++i){
    const double* ti = Tt.colptr(i);
    for(arma::uword k=0; k<nT; ++k){
      tSum(k) += ti[k];
    }
  }
  arma::uword nBlocks = countBlocks(nCol);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword colStart = blockStart(nCol, nBlocks, block);
    arma::uword colEnd = blockStart(nCol, nBlocks, block+1);
    // Four sets of 32 bins of nT values: in each, the low half byte's bins,
    // then the high's
    const arma::uword setLen = 32*nT;
    std::vector<double> bins(4*setLen), merged(setLen);
    for(arma::uword c=colStart; c<colEnd; ++c){
      const unsigned char* col = planes.colptr(c);
      std::fill(bins.begin(), bins.end(), 0.0);
      for(arma::uword i=0; i<nInd; ++i){
        const double* ti = Tt.colptr(i);
        // In a diploid each of an individual's two haplotypes has its own
        // set, and the next individual takes the other two
        double* pair = bins.data()+2*(i & 1u)*setLen;
        for(arma::uword h=0; h<ploidy; ++h){
          double* bh = (ploidy==2) ? pair+h*setLen : pair;
          unsigned int b = col[i*ploidy+h];
          double* lo = bh+(b & 15u)*nT;
          double* hi = bh+(16+(b>>4))*nT;
          for(arma::uword k=0; k<nT; ++k){
            lo[k] += ti[k];
            hi[k] += ti[k];
          }
        }
      }
      for(arma::uword e=0; e<setLen; ++e){
        merged[e] = (bins[e]+bins[setLen+e])+
          (bins[2*setLen+e]+bins[3*setLen+e]);
      }
      for(arma::uword b=0; b<8; ++b){
        arma::uword j = 8*c+b;
        if(j>=nLoci){
          break;
        }
        const double* mb = merged.data()+16*(b/4)*nT;
        unsigned int bit = 1u<<(b%4);
        double* wj = Wt.colptr(j);
        for(arma::uword k=0; k<nT; ++k){
          double sum = 0.0;
          for(unsigned int q=1; q<16; ++q){
            if(q & bit){
              sum += mb[q*nT+k];
            }
          }
          wj[k] = alleleScale*sum-offset(j)*tSum(k);
        }
      }
    }
  }
}

// Called by fastRRBLUP function for more than one trait
//
// The model is the multiple trait RR-BLUP model, with the marker effects of
// the traits correlated through Vu and the errors through Ve, and every trait
// recorded on every individual. The variance components are estimated by
// REML on a subset of the records, as for one trait, unless they are
// supplied.
//
// With every trait recorded and the same fixed effects for all of them, the
// canonical transformation separates the traits exactly. Writing Ve = L*L'
// and L^-1*Vu*L^-T = Q*diag(lam)*Q', the responses Y*T' for T = Q'*L^-1 have
// identity residual covariance and genetic covariance diag(lam), so trait k
// of the transformed responses is a single trait RR-BLUP model with variance
// ratio 1/lam_k. The marker effects are T^-1 times the transformed ones.
//
// The transformed models are solved together by preconditioned conjugate
// gradient, each trait by its own recurrence, with the products with the
// genotypes taken for every trait at once. The iteration works in whitened
// effects a_k = u_k/sqrt(lam_k), whose coefficient matrix is
// lam_k*Mc'*S*Mc+I. A direction with no genetic variance, lam_k of zero,
// then has no effects and needs no special case. The products, the
// absorption of the fixed effects and the preconditioner are kept apart from
// the iteration so that another residual structure only has to replace the
// last two; one that couples the traits would also need the recurrences
// joined.
//
// x is an indicator vector for the fixed effect levels, as in callRRBLUP
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP_MV(const arma::mat& Y, arma::uvec x,
                             const Rcpp::S4& pop,
                             arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                             const arma::mat& VuIn, const arma::mat& VeIn,
                             arma::uword maxIter, bool estVarComp,
                             arma::uvec subset, int nThreads){
  const arma::uword nT = Y.n_cols;
  if(!Y.is_finite()){
    Rcpp::stop("fastRRBLUP: every trait must be recorded on every individual");
  }
  PackedGeno g(pop);
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  arma::uword ploidy = g.ploidy;
  arma::uword nInd = g.nInd;
  arma::uword nLoci = pos.base.size();
  if(Y.n_rows!=nInd){
    Rcpp::stop("fastRRBLUP: the responses must have one row per individual");
  }
  arma::Mat<unsigned char> M = decodeChipDosage(g, pos, nThreads);
  double dP = double(ploidy);
  // Sum to zero design matrix with an intercept, matching callRRBLUP
  arma::mat Xfix = makeX(x);
  arma::uword q = Xfix.n_cols;
  arma::vec xa(ploidy+1);
  for(arma::uword i=0; i<xa.n_elem; ++i){
    xa(i) = (double(i)-dP/2.0)*(2.0/dP);
  }
  const double* xaPtr = xa.memptr();
  // The fixed effects are absorbed with S = I-X*(X'X)^-1*X', which is the
  // same for every trait because every trait has the same fixed effects
  arma::mat XtXi = inv_sympd(Xfix.t()*Xfix);
  arma::mat SY = Y-Xfix*(XtXi*(Xfix.t()*Y));
  arma::rowvec Xsum = sum(Xfix,0); // X'*1

  // Variance components, by REML on a subset of the records as for one
  // trait. Whichever of the two cross products is smaller is decomposed.
  arma::mat Vu, Ve;
  if(estVarComp){
    subset -= 1; // R to C++
    arma::uword nSub = subset.n_elem;
    arma::mat ySub = Y.rows(subset);
    arma::mat XSub = Xfix.rows(subset);
    arma::mat none;
    if(nLoci<nSub){
      arma::mat MSub(nSub,nLoci);
      for(arma::uword j=0; j<nLoci; ++j){
        const unsigned char* Mcol = M.colptr(j);
        double* MSubCol = MSub.colptr(j);
        for(arma::uword i=0; i<nSub; ++i){
          MSubCol[i] = xaPtr[Mcol[subset(i)]];
        }
      }
      remlMultiTrait(ySub, XSub, MSub, none, true, 1000, 1e-6, Vu, Ve);
    }else{
      arma::mat K = gramFromGeno(M, subset, xaPtr);
      remlMultiTrait(ySub, XSub, none, K, false, 1000, 1e-6, Vu, Ve);
    }
  }else{
    Vu = 0.5*(VuIn+VuIn.t());
    Ve = 0.5*(VeIn+VeIn.t());
  }
  if( (Vu.n_rows!=nT) || (Vu.n_cols!=nT) ||
      (Ve.n_rows!=nT) || (Ve.n_cols!=nT) ){
    Rcpp::stop("fastRRBLUP: Vu and Ve must be square with one row per trait");
  }

  // The canonical transformation
  arma::mat Le;
  if(!arma::chol(Le, Ve, "lower")){
    Rcpp::stop("fastRRBLUP: Ve must be positive definite");
  }
  arma::mat Linv;
  if(!arma::inv(Linv, arma::trimatl(Le))){
    Rcpp::stop("fastRRBLUP: Ve must be positive definite");
  }
  arma::mat A = Linv*Vu*Linv.t();
  A = 0.5*(A+A.t());
  arma::vec lam;
  arma::mat Qv;
  if(!arma::eig_sym(lam, Qv, A)){
    Rcpp::stop("fastRRBLUP: decomposition of Vu failed");
  }
  // These are the eigenvalues of a covariance matrix among traits, not of a
  // cross product of records, so rankCutoff does not apply. An eigenvalue
  // that is negative beyond rounding means Vu is not a covariance matrix, and
  // one that is negative within rounding is a direction with no genetic
  // variance. A positive eigenvalue is the variance ratio of its transformed
  // trait and is kept however small it is next to the largest, because the
  // ratios of different transformed traits need not be comparable: an error
  // variance close to zero for one trait makes its ratio huge without making
  // any other trait's ratio negligible.
  double lamMax = lam.max();
  if(!(lamMax>0.0)){
    Rcpp::stop("fastRRBLUP: estimated marker variance is not positive");
  }
  if(lam.min() < -1.0e-8*lamMax){
    Rcpp::stop("fastRRBLUP: Vu must be positive semidefinite");
  }
  lam.elem(find(lam<0.0)).zeros();
  arma::vec sqrtLam = sqrt(lam);
  arma::mat T = Qv.t()*Linv;
  arma::mat Tinv = Le*Qv;
  arma::mat SYs = SY*T.t(); // absorbed responses, transformed

  // One pass over the genotypes gives the column means, the right hand sides
  // and the diagonal of Mc'*S*Mc, from which every trait's preconditioner
  // follows. As for one trait, the pass splits over loci and accumulates
  // into scalars.
  arma::rowvec Mmean(nLoci);
  arma::vec dg(nLoci);
  arma::mat Rt(nT,nLoci); // right hand sides, one column per marker
  arma::rowvec SYsSum = sum(SYs,0);
  const double* XfixPtr = Xfix.memptr();
  {
    arma::uword nBlocks = countBlocks(nLoci);
    int setupThreads = nThreads;
    if(nBlocks < static_cast<arma::uword>(setupThreads)){
      setupThreads = static_cast<int>(nBlocks);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(setupThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword lociStart = blockStart(nLoci, nBlocks, block);
      arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
      arma::vec XtM(q), sSy(nT);
      for(arma::uword j=lociStart; j<lociEnd; ++j){
        const unsigned char* Mcol = M.colptr(j);
        double s1 = 0.0, s2 = 0.0;
        XtM.zeros();
        sSy.zeros();
        for(arma::uword i=0; i<nInd; ++i){
          double value = xaPtr[Mcol[i]];
          s1 += value;
          s2 += value*value;
          for(arma::uword k=0; k<nT; ++k){
            sSy(k) += value*SYs(i,k);
          }
          for(arma::uword k=0; k<q; ++k){
            XtM(k) += XfixPtr[i+k*nInd]*value;
          }
        }
        double mean = s1/double(nInd);
        Mmean(j) = mean;
        for(arma::uword k=0; k<nT; ++k){
          Rt(k,j) = sqrtLam(k)*(sSy(k)-mean*SYsSum(k));
        }
        // Centering the column shifts X'*Mc and takes n*mean^2 off Mc'*Mc
        for(arma::uword k=0; k<q; ++k){
          XtM(k) -= mean*Xsum(k);
        }
        double diagonal = s2-s1*mean-as_scalar(XtM.t()*XtXi*XtM);
        if(diagonal<0.0){
          // Rounding only, the true value cannot be negative
          diagonal = 0.0;
        }
        dg(j) = diagonal;
      }
    }
  }
  arma::mat precond = lam*dg.t()+1.0;

  // As for one trait, the iteration works from the markers packed to one bit
  // per allele, built only once the dosages have been freed
  M.reset();
  arma::Mat<unsigned char> planes = packChipPlanes(g, pos, nThreads);
  arma::vec offset = 1.0+Mmean.t();
  arma::mat tabs(32*nT,planes.n_cols);

  // Preconditioned conjugate gradient on the whitened effects. The
  // transformed traits are separate systems, so each has its own recurrence,
  // with its own step lengths and its own test of convergence against its
  // initial residual, as if it were fitted on its own. A shared recurrence
  // would let a trait whose residual is far larger, such as one with an
  // error variance close to zero, decide when every trait stops. The
  // recurrences share only the products with the genotypes, which are taken
  // for every trait at once. Each product with the coefficient matrix is
  //   sqrt(lam)*(Mc'*(S*(Mc*(sqrt(lam)*p)))) + p
  // with the scaling by sqrt(lam) applied trait by trait. A trait that has
  // converged keeps a search direction of zero, which leaves it unchanged.
  arma::mat At(nT,nLoci,arma::fill::zeros);
  arma::mat Zt = Rt/precond;
  arma::mat Pt = Zt;
  arma::mat Apt(nT,nLoci), Ps(nT,nLoci);
  arma::mat Tt(nT,nInd);
  arma::vec rNorm0 = sqrt(sum(square(Rt),1));
  arma::vec rz = sum(Rt%Zt,1);
  std::vector<bool> active(nT);
  arma::uword nActive = 0;
  for(arma::uword k=0; k<nT; ++k){
    active[k] = (rNorm0(k)>0.0);
    if(active[k]){
      ++nActive;
    }else{
      Pt.row(k).zeros();
    }
  }
  for(arma::uword iter=0; (iter<maxIter)&&(nActive>0); ++iter){
    Ps = Pt;
    Ps.each_col() %= sqrtLam;
    multMcBlock(Tt, planes, ploidy, offset, Ps, tabs, nThreads);
    Tt -= ((Tt*Xfix)*XtXi)*Xfix.t(); // S applied to each trait
    multMcTBlock(Apt, planes, ploidy, offset, Tt, nThreads);
    Apt.each_col() %= sqrtLam;
    Apt += Pt;
    arma::vec pAp = sum(Pt%Apt,1);
    for(arma::uword k=0; k<nT; ++k){
      if(!active[k]) continue;
      if(pAp(k)<=0.0){
        // The matrix is positive definite, so this is rounding once there is
        // nothing left to solve for
        active[k] = false;
        --nActive;
        Pt.row(k).zeros();
        continue;
      }
      double step = rz(k)/pAp(k);
      At.row(k) += step*Pt.row(k);
      Rt.row(k) -= step*Apt.row(k);
      if(norm(Rt.row(k),2)<=1.0e-8*rNorm0(k)){
        active[k] = false;
        --nActive;
        Pt.row(k).zeros();
        continue;
      }
      Zt.row(k) = Rt.row(k)/precond.row(k);
      double rzNew = accu(Rt.row(k)%Zt.row(k));
      Pt.row(k) = Zt.row(k)+(rzNew/rz(k))*Pt.row(k);
      rz(k) = rzNew;
    }
  }
  if(nActive>0){
    Rf_warning("callFastRRBLUP: reached maxIter without converging");
  }

  // Back to the marker effects of the original traits, U = U*T^-T, held here
  // with one column per marker
  At.each_col() %= sqrtLam;
  arma::mat Ut = Tinv*At;

  // Fixed effects, recovered from the solution. With the same fixed effects
  // for every trait and every trait recorded, generalized least squares
  // reduces to least squares trait by trait. They are found against the
  // centered genotypes, while the genetic values are predicted from the
  // uncentered ones, so the intercept takes off Mmean*U.
  multMcBlock(Tt, planes, ploidy, offset, Ut, tabs, nThreads);
  arma::mat B = XtXi*(Xfix.t()*(Y-Tt.t()));
  arma::mat alpha = Ut.t();
  arma::rowvec shift = Mmean*alpha;

  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")=-shift,
                            Rcpp::Named("mu")=B.row(0)-shift,
                            Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve);
}

// Called by RRBLUP function
// The markers of this function and of the other callRRBLUP functions below
// are read straight from the population's packed genotypes, see PackedGeno
// [[Rcpp::export]]
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x,
                      const Rcpp::S4& pop,
                      arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                      int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,0,ploidy,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  Rcpp::List ans = solveRRBLUP(y, X, M);
  arma::vec u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-as_scalar(Mmean*u),
                            Rcpp::Named("mu")=beta(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP2(arma::mat y, arma::uvec x,
                       const Rcpp::S4& pop,
                       arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                       double Vu, double Ve, double tol, int maxIter,
                       bool useEM, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,0,ploidy,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  Rcpp::List ans = solveRRBLUP_EM(y, X, M, Vu, Ve, 
                                  tol, maxIter, useEM);
  arma::vec u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-as_scalar(Mmean*u),
                            Rcpp::Named("mu")=beta(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Intercepts -Mmean*U, one per column of U. Each is taken as the product of a
// row with a column, which is how a single trait fit has always worked it
// out, so a single trait result is the same to the last bit.
arma::rowvec negMeanProduct(const arma::rowvec& Mmean, const arma::mat& U){
  arma::rowvec out(U.n_cols);
  for(arma::uword k=0; k<U.n_cols; ++k){
    out(k) = -as_scalar(Mmean*U.col(k));
  }
  return out;
}

// Called by RRBLUP_D function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_D(arma::mat y, arma::uvec x,
                        const Rcpp::S4& pop,
                        arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                        int maxIter, int nThreads){
  // Fit GS model
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::Mat<unsigned char> M = genoDosage(g,lociPerChr,lociLoc,0,ploidy,nThreads);
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Mlist(0));
  Mlist(1) = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Mlist(1),1));

  //Scaled genotypes
  double dP = double(ploidy);
  arma::vec xx(ploidy+1);
  for(arma::uword i=0; i<xx.n_elem; ++i)
    xx(i) = double(i);
  arma::vec xa = (xx-dP/2.0)*(2.0/dP);
  arma::vec xd = xx%(dP-xx)*(2.0/dP)*(2.0/dP);
  
  // Calculate genotype frequencies
  arma::mat freqMat(ploidy+1,M.n_cols);
  arma::uvec fixed(M.n_cols);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
    arma::uvec freq(ploidy+1,arma::fill::zeros);
    for(arma::uword i=0; i<M.n_rows; ++i){
      freq(M(i,j)) += 1;
    }
    if(any(freq == M.n_rows)){
      fixed(j) = 1;
    }else{
      fixed(j) = 0;
      for(arma::uword i=0; i<freq.n_elem; ++i){
        freqMat(i,j) = double(freq(i))/double(M.n_rows);
      }
    }
  }

  // The genotypes are only needed for the frequencies above, so they are
  // freed before the model is fitted rather than held through it
  arma::uword m = M.n_cols;
  M.reset();
  Rcpp::List ans = solveRRBLUPMK(y, X, Mlist, maxIter);

  // Clear memory
  y.reset();
  X.reset();
  Mlist.reset();

  // Solve for average effect, trait by trait. The frequencies above are the
  // same for every trait, and a multiple trait fit has one column of fixed
  // effects and of marker effects per trait.
  arma::mat beta = ans["beta"];
  arma::uword nT = beta.n_cols;
  arma::mat alpha(m,nT), d(m,nT);
  arma::field<arma::mat> u = ans["u"];
  for(arma::uword k=0; k<nT; ++k){
    double meanD = beta(beta.n_rows-1,k);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword j=0; j<m; ++j){
      double genoMu, gvMu;
      arma::vec gv(ploidy+1);
      if(fixed(j) == 1){
        d(j,k) = 0;
        alpha(j,k) = 0;
      }else{
        d(j,k) = u(1)(j,k) + meanD;
        gv = xa*u(0)(j,k) + xd*d(j,k);
        genoMu = accu(freqMat.col(j)%xa);
        gvMu = accu(freqMat.col(j)%gv);
        alpha(j,k) = accu(freqMat.col(j)%(gv-gvMu)%(xa-genoMu))/
          accu(freqMat.col(j)%(xa-genoMu)%(xa-genoMu));
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha")=alpha,
    Rcpp::Named("beta")=negMeanProduct(Mmean,alpha),
    Rcpp::Named("a")=u(0),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta.row(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}

// Called by RRBLUP_D2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_D2(arma::mat y, arma::uvec x,
                         const Rcpp::S4& pop,
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                         int maxIter, double Va, double Vd, double Ve,
                         double tol, bool useEM, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::Mat<unsigned char> M = genoDosage(g,lociPerChr,lociLoc,0,ploidy,nThreads);
  arma::mat Ma = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Ma);
  arma::mat Md = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Md,1));

  //Scaled genotypes
  double dP = double(ploidy);
  arma::vec xx(ploidy+1);
  for(arma::uword i=0; i<xx.n_elem; ++i)
    xx(i) = double(i);
  arma::vec xa = (xx-dP/2.0)*(2.0/dP);
  arma::vec xd = xx%(dP-xx)*(2.0/dP)*(2.0/dP);
  
  // Calculate genotype frequencies
  arma::mat freqMat(ploidy+1,M.n_cols);
  arma::uvec fixed(M.n_cols);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
    arma::uvec freq(ploidy+1,arma::fill::zeros);
    for(arma::uword i=0; i<M.n_rows; ++i){
      freq(M(i,j)) += 1;
    }
    if(any(freq == M.n_rows)){
      fixed(j) = 1;
    }else{
      fixed(j) = 0;
      for(arma::uword i=0; i<freq.n_elem; ++i){
        freqMat(i,j) = double(freq(i))/double(M.n_rows);
      }
    }
  }

  // The genotypes are only needed for the frequencies above, so they are
  // freed before the model is fitted rather than held through it
  arma::uword m = M.n_cols;
  M.reset();
  Rcpp::List ans = solveRRBLUP_EM2(y,X,Ma,Md,Va,Vd,Ve,tol,maxIter,useEM);

  // Clear memory
  y.reset();
  X.reset();
  Ma.reset();
  Md.reset();

  // Solve for average effect
  arma::vec alpha(m), d(m);
  arma::mat beta = ans["beta"];
  double meanD = beta(beta.n_elem-1);
  arma::mat u = ans["u"];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<m; ++j){
    double genoMu, gvMu;
    arma::vec gv(ploidy+1);
    if(fixed(j) == 1){
      d(j) = 0;
      alpha(j) = 0;
    }else{
      d(j) = u(j,1) + meanD;
      gv = xa*u(j,0) + xd*d(j);
      genoMu = accu(freqMat.col(j)%xa);
      gvMu = accu(freqMat.col(j)%gv);
      alpha(j) = accu(freqMat.col(j)%(gv-gvMu)%(xa-genoMu))/
        accu(freqMat.col(j)%(xa-genoMu)%(xa-genoMu));
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha")=alpha,
    Rcpp::Named("beta")=-as_scalar(Mmean*alpha),
    Rcpp::Named("a")=u.col(0),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_MV(arma::mat Y, arma::uvec x,
                         const Rcpp::S4& pop,
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                         int maxIter, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,0,ploidy,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  Rcpp::List ans = solveRRBLUPMV(Y, X, M, maxIter);
  arma::mat u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-(Mmean*u),
                            Rcpp::Named("mu")=beta.row(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Called by RRBLUP_GCA function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_GCA(arma::mat y, arma::uvec x,
                          const Rcpp::S4& pop,
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                          int maxIter, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,0,ploidy/2,
                         nThreads),ploidy/2,nThreads);
  Mlist(1) = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,ploidy/2,ploidy,
                         nThreads),ploidy/2,nThreads);
  arma::rowvec Mmean1 = mean(Mlist(0));
  arma::rowvec Mmean2 = mean(Mlist(1));
  arma::mat X = makeX(x);
  Rcpp::List ans = solveRRBLUPMK(y,X,Mlist,maxIter);
  // A multiple trait fit has one column of effects per trait
  arma::field<arma::mat> u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha1")=u(0),
                            Rcpp::Named("alpha2")=u(1),
                            Rcpp::Named("beta1")=negMeanProduct(Mmean1,u(0)),
                            Rcpp::Named("beta2")=negMeanProduct(Mmean2,u(1)),
                            Rcpp::Named("mu")=beta.row(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Called by RRBLUP_GCA2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_GCA2(arma::mat y, arma::uvec x,
                           const Rcpp::S4& pop,
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                           int maxIter, double Vu1, double Vu2, double Ve,
                           double tol, bool useEM, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::mat M1 = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,0,ploidy/2,
                                        nThreads),ploidy/2,nThreads);
  arma::mat M2 = genoToGenoA(genoDosage(g,lociPerChr,lociLoc,ploidy/2,ploidy,
                                        nThreads),ploidy/2,nThreads);
  arma::rowvec Mmean1 = mean(M1);
  arma::rowvec Mmean2 = mean(M2);
  arma::mat X = makeX(x);
  Rcpp::List ans = solveRRBLUP_EM2(y,X,M1,M2,Vu1,Vu2,Ve,tol,maxIter,useEM);
  arma::mat u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha1")=u.col(0),
                            Rcpp::Named("alpha2")=u.col(1),
                            Rcpp::Named("beta1")=-as_scalar(Mmean1*u.col(0)),
                            Rcpp::Named("beta2")=-as_scalar(Mmean2*u.col(1)),
                            Rcpp::Named("mu")=beta(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Called by RRBLUP_SCA function
// Only works with diploids
// [[Rcpp::export]]
Rcpp::List callRRBLUP_SCA(arma::mat y, arma::uvec x,
                          const Rcpp::S4& pop,
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                          int maxIter, int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  arma::field<arma::mat> Mlist(3);
  // The maternal and paternal dosages are each read once. Their sum is the
  // genotype dosage, the same byte for byte as reading it a third time.
  arma::Mat<unsigned char> M = genoDosage(g,lociPerChr,lociLoc,0,ploidy/2,
                                          nThreads);
  arma::Mat<unsigned char> Mpat = genoDosage(g,lociPerChr,lociLoc,ploidy/2,
                                             ploidy,nThreads);
  Mlist(0) = genoToGenoA(M,ploidy/2,nThreads);
  Mlist(1) = genoToGenoA(Mpat,ploidy/2,nThreads);
  M += Mpat;
  Mpat.reset();
  // Identify fixed markers
  arma::uword m = M.n_cols;
  double n = double(y.n_rows);
  arma::uvec fixed(m,arma::fill::ones);
  arma::rowvec p12(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
    p12(j) = accu(Mlist(0).col(j)%Mlist(1).col(j)+
      Mlist(0).col(j)+Mlist(1).col(j)+1)/(4*n);
    unsigned char firstGeno = M(0,j);
    for(arma::uword i=1; i<M.n_rows; ++i){
      if(firstGeno != M(i,j)){
        fixed(j) = 0;
        break;
      }
    }
  }
  
  // Set dominance genotypes
  Mlist(2) = genoToGenoD(M,ploidy,nThreads);
  M.reset();
  
  arma::rowvec Mmean1 = mean(Mlist(0));
  arma::rowvec Mmean2 = mean(Mlist(1));
  arma::mat X;
  X = join_rows(makeX(x),sum(Mlist(2),1));
  
  Rcpp::List ans = solveRRBLUPMK(y, X, Mlist, maxIter);
  
  // Clear memory
  y.reset();
  X.reset();
  Mlist.reset();
  
  // Solve for average effect, trait by trait. The frequencies above are the
  // same for every trait, and a multiple trait fit has one column of fixed
  // effects and of marker effects per trait.

  arma::rowvec p1 = (Mmean1+1)/2;
  arma::rowvec p2 = (Mmean2+1)/2;
  arma::rowvec intraD = p12-p1%p2;
  arma::mat beta = ans["beta"];
  arma::uword nT = beta.n_cols;
  arma::mat alpha1(m,nT), alpha2(m,nT), d(m,nT);
  arma::field<arma::mat> u = ans["u"];
  for(arma::uword k=0; k<nT; ++k){
    double meanD = beta(beta.n_rows-1,k);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword i=0; i<Mmean1.n_cols; ++i){
      double gv1, gv0, gvMu, alpha, numer, denom;
      if(fixed(i)==0){
        d(i,k) = u(2)(i,k) + meanD;
        //alpha1
        gv1 = (p2(i)+intraD(i)/p1(i))*(u(0)(i,k)+u(1)(i,k)) +
          ((1-p2(i))-intraD(i)/p1(i))*d(i,k);
        gv0 = (p2(i)-intraD(i)/(1-p1(i)))*d(i,k) +
          ((1-p2(i))+intraD(i)/(1-p1(i)))*(-u(0)(i,k)-u(1)(i,k));
        gvMu = p1(i)*gv1+(1-p1(i))*gv0;
        numer = p1(i)*(gv1-gvMu)*(1.0-Mmean1(i)) +
          (1-p1(i))*(gv0-gvMu)*(-1.0-Mmean1(i));
        denom = p1(i)*(1.0-Mmean1(i))*(1.0-Mmean1(i)) +
          (1-p1(i))*(-1.0-Mmean1(i))*(-1.0-Mmean1(i));
        alpha = numer/denom;
        if(!std::isfinite(alpha)) alpha=0;
        alpha1(i,k) = alpha;
        //alpha2
        gv1 = (p1(i)+intraD(i)/p2(i))*(u(0)(i,k)+u(1)(i,k)) +
          ((1-p1(i))-intraD(i)/p2(i))*d(i,k);
        gv0 = (p1(i)-intraD(i)/(1-p2(i)))*d(i,k) +
          ((1-p1(i))+intraD(i)/(1-p2(i)))*(-u(0)(i,k)-u(1)(i,k));
        gvMu = p2(i)*gv1+(1-p2(i))*gv0;
        numer = p2(i)*(gv1-gvMu)*(1.0-Mmean2(i)) +
          (1-p2(i))*(gv0-gvMu)*(-1.0-Mmean2(i));
        denom = p2(i)*(1.0-Mmean2(i))*(1.0-Mmean2(i)) +
          (1-p2(i))*(-1.0-Mmean2(i))*(-1.0-Mmean2(i));
        alpha = numer/denom;
        if(!std::isfinite(alpha)) alpha=0;
        alpha2(i,k) = alpha;
      }else{
        alpha1(i,k) = 0;
        alpha2(i,k) = 0;
        d(i,k) = 0;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha1")=alpha1,
    Rcpp::Named("alpha2")=alpha2,
    Rcpp::Named("beta1")=negMeanProduct(Mmean1,alpha1),
    Rcpp::Named("beta2")=negMeanProduct(Mmean2,alpha2),
    Rcpp::Named("a1")=u(0),
    Rcpp::Named("a2")=u(1),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta.row(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}

// Called by RRBLUP_SCA2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_SCA2(arma::mat y, arma::uvec x,
                           const Rcpp::S4& pop,
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                           int maxIter, double Vu1, double Vu2, double Vu3,
                           double Ve, double tol, bool useEM,
                           int nThreads){
  PackedGeno g(pop);
  arma::uword ploidy = g.ploidy;
  // The maternal and paternal dosages are each read once. Their sum is the
  // genotype dosage, the same byte for byte as reading it a third time.
  arma::Mat<unsigned char> M = genoDosage(g,lociPerChr,lociLoc,0,ploidy/2,
                                          nThreads);
  arma::Mat<unsigned char> Mpat = genoDosage(g,lociPerChr,lociLoc,ploidy/2,
                                             ploidy,nThreads);
  arma::mat M1 = genoToGenoA(M,ploidy/2,nThreads);
  arma::mat M2 = genoToGenoA(Mpat,ploidy/2,nThreads);
  M += Mpat;
  Mpat.reset();
  // Identify fixed markers
  arma::uword m = M.n_cols;
  double n = double(y.n_rows);
  arma::uvec fixed(m,arma::fill::ones);
  arma::rowvec p12(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
    p12(j) = accu(M1.col(j)%M2.col(j)+
      M1.col(j)+M2.col(j)+1)/(4*n);
    unsigned char firstGeno = M(0,j);
    for(arma::uword i=1; i<M.n_rows; ++i){
      if(firstGeno != M(i,j)){
        fixed(j) = 0;
        break;
      }
    }
  }
  
  // Set dominance genotypes
  arma::mat M3 = genoToGenoD(M,ploidy,nThreads);
  M.reset();
  
  arma::rowvec Mmean1 = mean(M1);
  arma::rowvec Mmean2 = mean(M2);
  arma::mat X;
  X = join_rows(makeX(x),sum(M3,1));
  
  Rcpp::List ans = solveRRBLUP_EM3(y,X,M1,M2,M3,Vu1,Vu2,Vu3,Ve,tol,maxIter,useEM);
  
  // Clear memory
  y.reset();
  X.reset();
  M1.reset();
  M2.reset();
  M3.reset();
  
  // Solve for average effect
  
  arma::rowvec p1 = (Mmean1+1)/2;
  arma::rowvec p2 = (Mmean2+1)/2;
  arma::rowvec intraD = p12-p1%p2;
  arma::vec alpha1(m), alpha2(m), d(m);
  arma::mat beta = ans["beta"];
  double meanD = beta(beta.n_elem-1);
  arma::mat u = ans["u"];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<Mmean1.n_cols; ++i){
    double gv1, gv0, gvMu, alpha, numer, denom;
    if(fixed(i)==0){
      d(i) = u(i,2) + meanD;
      //alpha1
      gv1 = (p2(i)+intraD(i)/p1(i))*(u(i,0)+u(i,1)) + 
        ((1-p2(i))-intraD(i)/p1(i))*d(i);
      gv0 = (p2(i)-intraD(i)/(1-p1(i)))*d(i) + 
        ((1-p2(i))+intraD(i)/(1-p1(i)))*(-u(i,0)-u(i,1));
      gvMu = p1(i)*gv1+(1-p1(i))*gv0;
      numer = p1(i)*(gv1-gvMu)*(1.0-Mmean1(i)) + 
        (1-p1(i))*(gv0-gvMu)*(-1.0-Mmean1(i));
      denom = p1(i)*(1.0-Mmean1(i))*(1.0-Mmean1(i)) + 
        (1-p1(i))*(-1.0-Mmean1(i))*(-1.0-Mmean1(i));
      alpha = numer/denom;
      if(!std::isfinite(alpha)) alpha=0;
      alpha1(i) = alpha;
      //alpha2
      gv1 = (p1(i)+intraD(i)/p2(i))*(u(i,0)+u(i,1)) + 
        ((1-p1(i))-intraD(i)/p2(i))*d(i);
      gv0 = (p1(i)-intraD(i)/(1-p2(i)))*d(i) + 
        ((1-p1(i))+intraD(i)/(1-p2(i)))*(-u(i,0)-u(i,1));
      gvMu = p2(i)*gv1+(1-p2(i))*gv0;
      numer = p2(i)*(gv1-gvMu)*(1.0-Mmean2(i)) + 
        (1-p2(i))*(gv0-gvMu)*(-1.0-Mmean2(i));
      denom = p2(i)*(1.0-Mmean2(i))*(1.0-Mmean2(i)) + 
        (1-p2(i))*(-1.0-Mmean2(i))*(-1.0-Mmean2(i));
      alpha = numer/denom;
      if(!std::isfinite(alpha)) alpha=0;
      alpha2(i) = alpha;
    }else{
      alpha1(i) = 0;
      alpha2(i) = 0;
      d(i) = 0;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha1")=alpha1,
    Rcpp::Named("alpha2")=alpha2,
    Rcpp::Named("beta1")=-as_scalar(Mmean1*alpha1),
    Rcpp::Named("beta2")=-as_scalar(Mmean2*alpha2),
    Rcpp::Named("a1")=u.col(0),
    Rcpp::Named("a2")=u.col(1),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}
