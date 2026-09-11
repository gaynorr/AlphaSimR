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

const double pi = 3.14159265358979323846;

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
class CoefMatFactor {
public:
  void compute(const arma::mat& C){
    arma::mat L;
    useChol = arma::chol(L, C, "lower");
    if(useChol){
      useChol = arma::inv(Linv, arma::trimatl(L));
    }
    if(useChol){
      Full.reset();
    }else{
      Linv.reset();
      Full = arma::inv(C);
    }
  }

  arma::mat solve(const arma::mat& b) const {
    if(useChol){
      return Linv.t()*(Linv*b);
    }
    return Full*b;
  }

  // Trace of the diagonal block of C^-1 spanning columns [first,last]
  double blockTrace(arma::uword first, arma::uword last) const {
    if(useChol){
      return arma::accu(arma::square(Linv.cols(first,last)));
    }
    return arma::sum(Full(arma::span(first,last),
                          arma::span(first,last)).diag());
  }

private:
  arma::mat Linv; // inverse of the lower Cholesky factor
  arma::mat Full; // full inverse, used only for the fallback
  bool useChol = false;
};

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

// Replacement for Armadillo's eig_sym
// Fixes an error with decomposition of large matrices
// If calcVec = false, eigvec is not used
// It would be better to template this function in the future
int eigen2(arma::vec& eigval, arma::mat& eigvec, arma::mat X,
           bool calcVec = true){
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

// Objective function for REML using the EMMA algorithm
Rcpp::List objREML(double param, Rcpp::List args){
  double df = args["df"];
  arma::vec eta = args["eta"];
  arma::vec lambda = args["lambda"];
  double value = df * log(sum(eta%eta/(lambda+param)));
  value += sum(log(lambda+param));
  return Rcpp::List::create(Rcpp::Named("objective") = value,
                            Rcpp::Named("output") = 0);
}

// Objective function for REML using the EMMA algorithm, for a system whose
// random effect is rank deficient. The eigenvalues that are zero are handled
// analytically rather than being stored: nNull of them contribute
// nNull*log(param) to the log determinant and Rnull/param to the quadratic
// form. Setting Rnull and nNull to zero recovers objREML.
Rcpp::List objREML2(double param, Rcpp::List args){
  double df = args["df"];
  arma::vec eta = args["eta"];
  arma::vec lambda = args["lambda"];
  double Rnull = args["Rnull"];
  double nNull = args["nNull"];
  double value = df * log(sum(eta%eta/(lambda+param)) + Rnull/param);
  value += sum(log(lambda+param)) + nNull*log(param);
  return Rcpp::List::create(Rcpp::Named("objective") = value,
                            Rcpp::Named("output") = 0);
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
  const int nGrid = 100;
  double logLower = log(1.0e-10);
  double logUpper = log(1.0e10);
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
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e} using the
//' EMMA algorithm \insertCite{kang_2008}{AlphaSimR}.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M){
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  arma::uword m = M.n_cols;
  double df = double(n)-double(q);

  arma::mat XtXi = inv_sympd(X.t()*X);
  arma::vec yv = y.col(0);

  // Response with the fixed effects absorbed, S*y for the projector
  // S = I-X*(X'X)^-1*X'. S is never formed.
  arma::vec Sy = yv - X*(XtXi*(X.t()*yv));
  double ySy = accu(Sy%Sy);

  // REML needs the nonzero eigenvalues of S*M*M'*S and the projections of
  // S*y onto their eigenvectors. Note that S*M*M'*S = A*A' for A = S*M, and
  // that A costs only O(n*q*m) because S is a projector, so neither M*M' nor
  // a multiplication by S is needed to reach them.
  arma::vec lambda; // nonzero eigenvalues
  arma::vec eta;    // matching projections of S*y
  arma::mat H;      // M*M', only formed when that branch needs it

  if(m<n){
    // A'*A and A*A' share their nonzero eigenvalues, so decompose the
    // smaller of the two
    arma::mat A = M - X*(XtXi*(X.t()*M));
    arma::vec eigval(m);
    arma::mat eigvec(m,m);
    eigen2(eigval, eigvec, A.t()*A);
    double maxEig = eigval.max();
    if(maxEig<0.0) maxEig = 0.0;
    arma::uvec keep = find(eigval>(double(m)*arma::datum::eps*maxEig));
    if(keep.n_elem>arma::uword(df)) keep = keep.tail(arma::uword(df));
    lambda = eigval(keep);
    // The eigenvectors of A*A' are A*w/sqrt(lambda), so the projections
    // follow from A'*S*y without ever forming them
    eta = (eigvec.cols(keep).t()*(A.t()*Sy))/sqrt(lambda);
  }else{
    H = M*M.t();
    // S*H*S expanded with P = X*(X'X)^-1*X'. Every term is O(q*n^2),
    // against O(n^3) for two dense multiplications by S.
    arma::mat G;
    {
      arma::mat PH = X*(XtXi*(X.t()*H));
      G = H - PH - PH.t() + ((PH*X)*XtXi)*X.t();
    }
    arma::vec eigval(n);
    arma::mat eigvec(n,n);
    eigen2(eigval, eigvec, G);
    G.reset();
    double maxEig = eigval.max();
    if(maxEig<0.0) maxEig = 0.0;
    arma::uvec keep = find(eigval>(double(n)*arma::datum::eps*maxEig));
    if(keep.n_elem>arma::uword(df)) keep = keep.tail(arma::uword(df));
    lambda = eigval(keep);
    eta = eigvec.cols(keep).t()*Sy;
  }

  // The remaining directions in the range of S have a zero eigenvalue. They
  // are summarized by their count and their total sum of squares.
  double Rnull = ySy-accu(eta%eta);
  if(Rnull<0.0) Rnull = 0.0;
  double nNull = df-double(lambda.n_elem);

  // Estimate variances
  Rcpp::List optRes = optimize(*objREML2,
                               Rcpp::List::create(
                                 Rcpp::Named("df")=df,
                                 Rcpp::Named("eta")=eta,
                                 Rcpp::Named("lambda")=lambda,
                                 Rcpp::Named("Rnull")=Rnull,
                                 Rcpp::Named("nNull")=nNull),
                                 1.0e-10, 1.0e10);
  double delta = optRes["parameter"];

  // Solve the mixed model equations. V = M*M'+delta*I is applied to X and y
  // without ever being inverted.
  arma::mat VX, Vy;
  if(m<n){
    // Woodbury identity, so the largest matrix factorized is m by m
    arma::mat C = M.t()*M;
    C.diag() += delta;
    arma::mat Ci = inv_sympd(C);
    VX = (X - M*(Ci*(M.t()*X)))/delta;
    Vy = (y - M*(Ci*(M.t()*y)))/delta;
  }else{
    H.diag() += delta;
    arma::mat sol;
    if(!solve(sol, H, join_rows(X,y))){
      Rcpp::stop("solveRRBLUP: mixed model equations could not be solved");
    }
    VX = sol.cols(0,q-1);
    Vy = sol.col(q);
  }
  arma::mat beta = solve(X.t()*VX, X.t()*Vy);
  arma::mat u = M.t()*(Vy-VX*beta);

  double Vu = (sum(eta%eta/(lambda+delta))+Rnull/delta)/df;
  double Ve = delta*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u);
}

//' @title Solve RR-BLUP with FaST-LMM
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e}. Takes the
//' same arguments and returns the same values as \code{\link{solveRRBLUP}},
//' but solves the mixed model equations using the factored spectral approach
//' of FaST-LMM \insertCite{lippert_2011}{AlphaSimR} rather than the EMMA
//' algorithm \insertCite{kang_2008}{AlphaSimR}. It is intended as an
//' eventual replacement for \code{\link{solveRRBLUP}}.
//'
//' @details
//' Both algorithms reduce the mixed model to a one dimensional search over
//' delta, the ratio of the residual variance to the marker variance. They
//' differ in the decomposition they search over.
//'
//' EMMA works with the fixed effects projected out. It needs the nonzero
//' eigenvalues of S*M*M'*S for the projector S, and then has to factorize
//' M*M'+delta*I a second time to reach the solutions.
//'
//' FaST-LMM works with M*M' itself and re-estimates the fixed effects at
//' every delta. One decomposition therefore supplies the likelihood, the
//' generalized least squares solution for the fixed effects, and the BLUPs,
//' and that decomposition is taken in whichever of the two spaces is
//' smaller. When there are fewer markers than records the eigenvectors of
//' M'*M serve in place of those of M*M', the rank deficient directions are
//' summarized analytically, and no matrix larger than M is ever formed.
//'
//' This is an independent implementation of the published method. It shares
//' no code with the FaST-LMM software distributed by Microsoft.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP2(const arma::mat& y, const arma::mat& X,
                        const arma::mat& M){
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  arma::uword m = M.n_cols;
  double df = double(n)-double(q);

  arma::vec yv = y.col(0);

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
    arma::mat MtX = M.t()*X;
    arma::vec eigval(m);
    arma::mat eigvec(m,m);
    eigen2(eigval, eigvec, M.t()*M);
    double maxEig = eigval.max();
    if(maxEig<0.0) maxEig = 0.0;
    arma::uvec keep = find(eigval>(double(m)*arma::datum::eps*maxEig));
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
    yy0 = accu(yv%yv)-accu(yt%yt);
    if(yy0<0.0) yy0 = 0.0;
    Xy0 = X.t()*yv-Xt.t()*yt;
    XX0 = X.t()*X-Xt.t()*Xt;
  }else{
    // M*M' is the smaller of the two, and its decomposition spans all of
    // R^n, so there is no null space to account for
    arma::mat K = M*M.t();
    arma::vec eigval(n);
    U.set_size(n,n);
    eigen2(eigval, U, K);
    K.reset();
    eigval.elem(find(eigval<0.0)).zeros(); // rounding only
    s2 = eigval;
    yt = U.t()*yv;
    Xt = U.t()*X;
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
    Rcpp::stop("solveRRBLUP2: mixed model equations could not be solved");
  }
  arma::mat beta = fit.beta;

  // The BLUPs are u = M'*H^-1*(y-X*beta), which the decomposition already
  // in hand applies without a second factorization
  arma::vec et = yt-Xt*fit.beta; // U'*(y-X*beta)
  arma::vec etd = et/(s2+delta);
  arma::mat u;
  if(m<n){
    arma::vec Mte = M.t()*(yv-X*fit.beta);
    u = MtU*etd+(Mte-MtU*et)/delta;
  }else{
    u = M.t()*(U*etd);
  }

  double Vu = fit.r/df;
  double Ve = delta*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u);
}

//' @title Solve Multivariate RR-BLUP
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Mu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X,
                         const arma::mat& M, int maxIter=1000, 
                         double tol=1e-6){
  arma::uword n = Y.n_rows;
  arma::uword m = Y.n_cols;
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, M*M.t());
  arma::mat Yt = Y.t()*eigvec;
  arma::mat Xt = X.t()*eigvec;
  arma::mat Vu = cov(Y)/2;
  arma::mat Ve = Vu;
  arma::mat W = Xt.t()*inv_sympd(Xt*Xt.t());
  arma::mat B = Yt*W; //BLUEs
  arma::cube Hinv(m,m,n);
  arma::mat tolEye = tol*arma::eye(m,m);
  arma::mat Gt(m,n), sigma(m,m), BNew, 
  VeNew(m,m), VuNew(m,m);
  double denom, numer;
  bool converging=true;
  int iter=0;
  while(converging){
    ++iter;
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    // The same m by m inverse is needed again in the second loop and
    // for the BLUPs, so it is formed once per iteration and kept
    for(arma::uword i=0; i<n; ++i){
      Hinv.slice(i) = inv_sympd(eigval(i)*Vu+Ve+tolEye);
      Gt.col(i) = (eigval(i)*Vu)*Hinv.slice(i)*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(arma::uword i=0; i<n; ++i){
      arma::mat lVu = eigval(i)*Vu;
      sigma = lVu-lVu*Hinv.slice(i)*lVu;
      arma::mat resid = Yt.col(i)-BNew*Xt.col(i)-Gt.col(i);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*(resid*resid.t()+sigma);
    }
    denom = fabs(sum(Ve.diag()));
    if(denom>0.0){
      numer = fabs(sum(VeNew.diag()-Ve.diag()));
      if((numer/denom)<tol) converging=false;
    }
    Ve = VeNew;
    Vu = VuNew;
    B = BNew;
    if(iter>=maxIter){
      Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
      break;
    }
  }
  // The BLUPs follow from the eigendecomposition already computed. The
  // variance of vec(E) is (eigvec (x) I)*blockdiag(eigval_i*Vu+Ve)*
  // (eigvec' (x) I), so its inverse is applied by rotating, solving each
  // m by m block, and rotating back. The n*m by n*m matrix that the
  // Kronecker form would build is never needed.
  for(arma::uword i=0; i<n; ++i){
    Hinv.slice(i) = inv_sympd(eigval(i)*Vu+Ve+tolEye);
  }
  arma::mat E = Y.t() - B*X.t();
  arma::mat Et = E*eigvec;
  arma::mat Wt(m,n);
  for(arma::uword i=0; i<n; ++i){
    Wt.col(i) = Hinv.slice(i)*Et.col(i);
  }
  arma::mat U = Vu*(Wt*eigvec.t())*M; //BLUPs
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multikernel RR-BLUP
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Mlist a list of M matrices
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMK(arma::mat& y, arma::mat& X,
                         arma::field<arma::mat>& Mlist,
                         int maxIter=40){
  double tol = 1e-4;
  double minVC = sqrt(2.2204460492503131e-016);
  arma::uword k = Mlist.n_elem;
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  double df = double(n)-double(q);
  arma::field<arma::mat> V(k);
  for(arma::uword i=0; i<k; ++i){
    V(i) = Mlist(i)*Mlist(i).t();
  }
  // Choose once between forming the n by n products WQX*V(i) and working
  // from the thin matrices in Mlist. Both give the same average information
  // matrix, but the second is cheaper only when the number of columns is
  // small relative to n, roughly below n/2, so the flop counts decide.
  double denseCost=0.0, thinCost=0.0;
  for(arma::uword i=0; i<k; ++i){
    double mi = double(Mlist(i).n_cols);
    denseCost += double(n)*double(n)*double(n);
    thinCost += double(n)*double(n)*mi;
    for(arma::uword j=i; j<k; ++j){
      denseCost += double(n)*double(n);
      thinCost += double(n)*mi*double(Mlist(j).n_cols);
    }
  }
  bool useThin = thinCost<denseCost;
  arma::field<arma::mat> C(k);
  arma::mat A(k+1,k+1), W0(n,n), W(n,n), WX(n,q), WQX(n,n);
  arma::vec qvec(k+1), sigma(k+1);
  arma::mat g(n,1);
  double logdetV=0;
  double rss, ldet, llik, llik0=0, deltaLlik, taper, 
  value, sign;
  bool invPass;
  arma::field<arma::mat> T(k);
  sigma.fill(var(y.col(0)));
  int iter = 0;
  while(true){
    ++iter;
    W0 = V(0)*sigma(0);
    for(arma::uword i=1; i<k; ++i){
      W0 += V(i)*sigma(i);
    }
    W0.diag() += sigma(k);
    invPass = invAndLogDet(W,logdetV,W0);
    if(!invPass){
      W = pinv(W0);
      log_det(value, sign, W0);
      logdetV = value*sign;
    }
    WX = W*X;
    WQX = W - WX*solve(X.t()*WX, WX.t());
    g = WQX*y;
    rss = as_scalar(y.t()*g);
    sigma = sigma*(rss/df);
    WQX = WQX*(df/rss);
    g = g*(df/rss);
    // WQX has rank n-q, so its ordinary determinant is zero and
    // log_det cannot be used on it. What the likelihood needs is the
    // product of its nonzero eigenvalues, which follows from
    // |WQX|+ = |X'X|/(|V|*|X'V^-1*X|). The constant |X'X| is dropped
    // because only changes in llik are used.
    log_det(value, sign, X.t()*WX);
    ldet = -logdetV - value*sign + df*log(df/rss);
    llik = ldet/2 - df/2;
    if(iter == 1) llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    // y'*T_i*WQX*y equals (WQX*y)'*V_i*(WQX*y), which is quadratic
    // rather than cubic in n. A is symmetric, so only its upper
    // triangle is computed, and WQX is symmetric, so the transposes
    // that would each need a temporary n by n matrix are dropped.
    if(useThin){
      // V(i) = Mlist(i)*Mlist(i)', so every trace the average information
      // matrix needs is reachable through C(i) = WQX*Mlist(i), without ever
      // forming the n by n product WQX*V(i):
      //   tr(WQX*V_i)         = accu(Mlist(i) % C(i))
      //   tr(WQX*V_i*WQX*V_j) = accu(B % B) for B = Mlist(i).t()*C(j)
      //   tr(WQX*V_i*WQX)     = accu(C(i) % C(i))
      for(arma::uword i=0; i<k; ++i){
        C(i) = WQX*Mlist(i);
      }
      for(arma::uword i=0; i<k; ++i){
        qvec(i) = accu(square(Mlist(i).t()*g)) - accu(Mlist(i)%C(i));
        for(arma::uword j=i; j<k; ++j){
          arma::mat B = Mlist(i).t()*C(j);
          A(i,j) = accu(B%B);
          A(j,i) = A(i,j);
        }
        A(i,k) = accu(C(i)%C(i));
        A(k,i) = A(i,k);
      }
    }else{
      for(arma::uword i=0; i<k; ++i){
        T(i) = WQX*V(i);
      }
      for(arma::uword i=0; i<k; ++i){
        qvec(i) = as_scalar(g.t()*V(i)*g) - sum(T(i).diag());
        for(arma::uword j=i; j<k; ++j){
          A(i,j) = accu(T(i)%T(j).t());
          A(j,i) = A(i,j);
        }
        A(i,k) = accu(T(i)%WQX);
        A(k,i) = A(i,k);
      }
    }
    A(k,k) = accu(WQX%WQX);
    qvec(k) = as_scalar(g.t()*g) - sum(WQX.diag());
    A = pinv(A);
    qvec = A*qvec;
    if(iter == 1){
      taper = 0.5;
    }else if(iter == 2){
      taper = 0.7;
    }else{
      taper = 0.9;
    }
    sigma += taper*qvec;
    while(sigma.min() < minVC){
      sigma(sigma.index_min()) = minVC;
    }
    if((iter>1)&(fabs(deltaLlik)<tol*10)){
      break;
    }
    if(max(abs(qvec)) < tol){
      break;
    }
    if(iter >= maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  arma::mat beta(q,1), ee(n,1);
  arma::field<arma::mat> u(k);
  beta = solve(X.t()*W*X,X.t()*W*y);
  ee = y - X*beta;
  for(arma::uword i=0; i<k; ++i){
    u(i) = sigma(i)*Mlist(i).t()*W*ee;
  }
  arma::vec Vu(k), Ve(1);
  Vu = sigma(arma::span(0,k-1));
  Ve = sigma(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("iter")=iter);
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
  if(!useEM & (n<m)){
    arma::mat Vinv = inv_sympd(M*M.t()*Vu+arma::eye(n,n)*Ve);
    arma::mat beta = solve(X.t()*Vinv*X, X.t()*Vinv*Y);
    arma::mat u = M.t()*Vinv*(Y-X*beta)*Vu;
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
  RHSfac.compute(RHS);
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
  RHSfac.compute(RHS);
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
//' @param Vu3 initial guess for variance of the second marker effects
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
  RHSfac.compute(RHS);
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

// Applies the centered genotype matrix to a vector, t = Mc*p.
//
// The dosages are never held anywhere. A column of M is one byte per record
// and is turned into a dosage by a lookup as it is read, which is what keeps
// the memory used by the solver close to the size of the genotypes.
//
// Centering is a rank one correction, Mc*p = M*p - 1*(mean'p), so no centered
// copy of the genotypes is needed either.
//
// The work is split over blocks of records. A block owns its own part of t,
// so nothing is accumulated across threads and the result does not depend on
// how many threads are used.
void multMc(arma::vec& t, const arma::Mat<unsigned char>& M,
            const double* xaPtr, const arma::rowvec& Mmean,
            const arma::vec& p, int nThreads){
  arma::uword nInd = M.n_rows;
  arma::uword nLoci = M.n_cols;
  const double* pPtr = p.memptr();
  double shift = 0.0;
  for(arma::uword j=0; j<nLoci; ++j){
    shift += Mmean(j)*pPtr[j];
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
    for(arma::uword j=0; j<nLoci; ++j){
      const unsigned char* Mcol = M.colptr(j);
      double pj = pPtr[j];
      for(arma::uword i=indStart; i<indEnd; ++i){
        tPtr[i] += xaPtr[Mcol[i]]*pj;
      }
    }
  }
}

// Applies the transpose of the centered genotype matrix, w = Mc'*t.
//
// Split over blocks of loci. A block works out its own entries of w from its
// own columns of M, so again nothing is accumulated across threads.
void multMcT(arma::vec& w, const arma::Mat<unsigned char>& M,
             const double* xaPtr, const arma::rowvec& Mmean,
             const arma::vec& t, int nThreads){
  arma::uword nInd = M.n_rows;
  arma::uword nLoci = M.n_cols;
  const double* tPtr = t.memptr();
  double tSum = 0.0;
  for(arma::uword i=0; i<nInd; ++i){
    tSum += tPtr[i];
  }
  double* wPtr = w.memptr();
  arma::uword nBlocks = countBlocks(nLoci);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword lociStart = blockStart(nLoci, nBlocks, block);
    arma::uword lociEnd = blockStart(nLoci, nBlocks, block+1);
    for(arma::uword j=lociStart; j<lociEnd; ++j){
      const unsigned char* Mcol = M.colptr(j);
      double dot = 0.0;
      for(arma::uword i=0; i<nInd; ++i){
        dot += xaPtr[Mcol[i]]*tPtr[i];
      }
      wPtr[j] = dot-Mmean(j)*tSum;
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
// using the same FaST-LMM machinery as solveRRBLUP2. Only the variance
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
  eigval.elem(find(eigval<0.0)).zeros(); // rounding only
  arma::vec yt = U.t()*yv;
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
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                          double Vu, double Ve, arma::uword maxIter,
                          bool estVarComp, arma::uvec subset, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::uword nInd = M.n_rows;
  arma::uword nLoci = M.n_cols;
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
  // there are fewer markers than sampled records, solveRRBLUP2 is given the
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
      Rcpp::List varComp = solveRRBLUP2(ySub, XSub, MSub);
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
    multMc(t, M, xaPtr, Mmean, pDir, nThreads);
    t -= Xfix*(XtXi*(Xfix.t()*t)); // S*t
    multMcT(Ap, M, xaPtr, Mmean, t, nThreads);
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
  
  // Fixed effects, recovered from the solution
  multMc(t, M, xaPtr, Mmean, u, nThreads);
  arma::vec beta = XtXi*(Xfix.t()*(y-t));
  
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-as_scalar(Mmean*u),
                            Rcpp::Named("mu")=beta(0),
                            Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve);
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x,
                      arma::field<arma::Cube<unsigned char> >& geno, 
                      arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                      int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
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
                       arma::field<arma::Cube<unsigned char> >& geno, 
                       arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                       double Vu, double Ve, double tol, int maxIter, 
                       bool useEM, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
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

// Called by RRBLUP_D function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_D(arma::mat y, arma::uvec x,
                        arma::field<arma::Cube<unsigned char> >& geno, 
                        arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                        int maxIter, int nThreads){
  // Fit GS model
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Mlist(0));
  Mlist(1) = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Mlist(1),1));
  Rcpp::List ans = solveRRBLUPMK(y, X, Mlist, maxIter);
  
  // Clear memory
  y.reset();
  X.reset();
  Mlist.reset();
  
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
  
  // Solve for average effect
  arma::vec alpha(M.n_cols), d(M.n_cols);
  arma::mat beta = ans["beta"];
  double meanD = beta(beta.n_elem-1);
  arma::field<arma::mat> u = ans["u"];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
    double genoMu, gvMu;
    arma::vec gv(ploidy+1);
    if(fixed(j) == 1){
      d(j) = 0;
      alpha(j) = 0;
    }else{
      d(j) = u(1).at(j) + meanD;
      gv = xa*u(0).at(j) + xd*d(j);
      genoMu = accu(freqMat.col(j)%xa);
      gvMu = accu(freqMat.col(j)%gv);
      alpha(j) = accu(freqMat.col(j)%(gv-gvMu)%(xa-genoMu))/
        accu(freqMat.col(j)%(xa-genoMu)%(xa-genoMu));
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha")=alpha,
    Rcpp::Named("beta")=-as_scalar(Mmean*alpha),
    Rcpp::Named("a")=u(0),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}

// Called by RRBLUP_D2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_D2(arma::mat y, arma::uvec x, 
                         arma::field<arma::Cube<unsigned char> >& geno, 
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                         int maxIter, double Va, double Vd, double Ve, 
                         double tol, bool useEM, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::mat Ma = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Ma);
  arma::mat Md = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Md,1));
  Rcpp::List ans = solveRRBLUP_EM2(y,X,Ma,Md,Va,Vd,Ve,tol,maxIter,useEM);
  
  // Clear memory
  y.reset();
  X.reset();
  Ma.reset();
  Md.reset();
  
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
    if(any(freq == M.n_cols)){
      fixed(j) = 1;
    }else{
      fixed(j) = 0;
      for(arma::uword i=0; i<freq.n_elem; ++i){
        freqMat(i,j) = double(freq(i))/double(M.n_rows);
      }
    }
  }
  
  // Solve for average effect
  arma::vec alpha(M.n_cols), d(M.n_cols);
  arma::mat beta = ans["beta"];
  double meanD = beta(beta.n_elem-1);
  arma::mat u = ans["u"];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword j=0; j<M.n_cols; ++j){
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
                         arma::field<arma::Cube<unsigned char> >& geno, 
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                         int maxIter, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
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
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                          int maxIter, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  Mlist(1) = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  arma::rowvec Mmean1 = mean(Mlist(0));
  arma::rowvec Mmean2 = mean(Mlist(1));
  arma::mat X = makeX(x);
  Rcpp::List ans = solveRRBLUPMK(y,X,Mlist,maxIter);
  arma::field<arma::mat> u = ans["u"];
  arma::mat beta = ans["beta"];
  return Rcpp::List::create(Rcpp::Named("alpha1")=u(0),
                            Rcpp::Named("alpha2")=u(1),
                            Rcpp::Named("beta1")=-as_scalar(Mmean1*u(0)),
                            Rcpp::Named("beta2")=-as_scalar(Mmean2*u(1)),
                            Rcpp::Named("mu")=beta(0),
                            Rcpp::Named("Vu")=ans["Vu"],
                            Rcpp::Named("Ve")=ans["Ve"]);
}

// Called by RRBLUP_GCA2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_GCA2(arma::mat y, arma::uvec x, 
                           arma::field<arma::Cube<unsigned char> >& geno, 
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                           int maxIter, double Vu1, double Vu2, double Ve, 
                           double tol, bool useEM, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat M1 = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                                             nThreads),ploidy/2,nThreads);
  arma::mat M2 = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
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
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                          int maxIter, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::field<arma::mat> Mlist(3);
  Mlist(0) = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  Mlist(1) = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  // Identify fixed markers
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,
                                       lociLoc, nThreads);
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
  arma::rowvec Mmean12 = mean(Mlist(0)%Mlist(1));
  arma::mat X;
  X = join_rows(makeX(x),sum(Mlist(2),1));
  
  Rcpp::List ans = solveRRBLUPMK(y, X, Mlist, maxIter);
  
  // Clear memory
  y.reset();
  X.reset();
  Mlist.reset();
  
  // Solve for average effect
  
  arma::rowvec p1 = (Mmean1+1)/2;
  arma::rowvec p2 = (Mmean2+1)/2;
  arma::rowvec intraD = p12-p1%p2;
  arma::vec alpha1(m), alpha2(m), d(m);
  arma::mat beta = ans["beta"];
  double meanD = beta(beta.n_elem-1);
  arma::field<arma::mat> u = ans["u"];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword i=0; i<Mmean1.n_cols; ++i){
    double gv1, gv0, gvMu, alpha, numer, denom;
    if(fixed(i)==0){
      d(i) = u(2).at(i) + meanD;
      //alpha1
      gv1 = (p2(i)+intraD(i)/p1(i))*(u(0).at(i)+u(1).at(i)) + 
        ((1-p2(i))-intraD(i)/p1(i))*d(i);
      gv0 = (p2(i)-intraD(i)/(1-p1(i)))*d(i) + 
        ((1-p2(i))+intraD(i)/(1-p1(i)))*(-u(0).at(i)-u(1).at(i));
      gvMu = p1(i)*gv1+(1-p1(i))*gv0;
      numer = p1(i)*(gv1-gvMu)*(1.0-Mmean1(i)) + 
        (1-p1(i))*(gv0-gvMu)*(-1.0-Mmean1(i));
      denom = p1(i)*(1.0-Mmean1(i))*(1.0-Mmean1(i)) + 
        (1-p1(i))*(-1.0-Mmean1(i))*(-1.0-Mmean1(i));
      alpha = numer/denom;
      if(!std::isfinite(alpha)) alpha=0;
      alpha1(i) = alpha;
      //alpha2
      gv1 = (p1(i)+intraD(i)/p2(i))*(u(0).at(i)+u(1).at(i)) + 
        ((1-p1(i))-intraD(i)/p2(i))*d(i);
      gv0 = (p1(i)-intraD(i)/(1-p2(i)))*d(i) + 
        ((1-p1(i))+intraD(i)/(1-p2(i)))*(-u(0).at(i)-u(1).at(i));
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
    Rcpp::Named("a1")=u(0),
    Rcpp::Named("a2")=u(1),
    Rcpp::Named("d")=d,
    Rcpp::Named("mu")=beta(0),
    Rcpp::Named("Vu")=ans["Vu"],
    Rcpp::Named("Ve")=ans["Ve"]
  );
}

// Called by RRBLUP_SCA2 function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_SCA2(arma::mat y, arma::uvec x, 
                           arma::field<arma::Cube<unsigned char> >& geno, 
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                           int maxIter, double Vu1, double Vu2, double Vu3, 
                           double Ve, double tol, bool useEM,  
                           int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat M1 = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                                             nThreads),ploidy/2,nThreads);
  arma::mat M2 = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
                                             nThreads),ploidy/2,nThreads);
  // Identify fixed markers
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,
                                       lociLoc, nThreads);
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
  arma::rowvec Mmean12 = mean(M1%M2);
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

//' @title Solve Univariate Model
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Zu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K){
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));
  
  // Construct system of equations for eigendecomposition
  arma::mat XtXi = inv_sympd(X.t()*X);
  arma::mat ZK = Z*K;
  arma::mat H = ZK*Z.t(); // Used later
  H.diag() += offset;
  // S*H*S for the projector S = I-X*(X'X)^-1*X', expanded in terms of
  // P = X*(X'X)^-1*X'. Every term is O(q*n^2), against O(n^3) for two dense
  // multiplications by S, and S itself never has to be formed.
  arma::mat S;
  {
    arma::mat PH = X*(XtXi*(X.t()*H));
    S = H - PH - PH.t() + ((PH*X)*XtXi)*X.t();
  }
  
  // Compute eigendecomposition
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, S);
  S.reset();
  
  // Drop eigenvalues
  eigval = eigval(arma::span(q,eigvec.n_cols-1)) - offset;
  eigvec = eigvec(arma::span(0,eigvec.n_rows-1),
                  arma::span(q,eigvec.n_cols-1));
  
  // Estimate variances and solve equations
  arma::vec eta = eigvec.t()*y;
  Rcpp::List optRes = optimize(*objREML,
                               Rcpp::List::create(
                                 Rcpp::Named("df")=df,
                                 Rcpp::Named("eta")=eta,
                                 Rcpp::Named("lambda")=eigval),
                                 1.0e-10, 1.0e10);
  double delta = optRes["parameter"];
  H.diag() += (delta-offset);
  // V = Z*K*Z'+delta*I is applied to X and y rather than inverted. One
  // factorization serves both, and V^-1*(y-X*beta) = Vy-VX*beta follows by
  // linearity, so no second solve is needed.
  arma::mat sol;
  if(!solve(sol, H, join_rows(X,y))){
    Rcpp::stop("solveUVM: mixed model equations could not be solved");
  }
  arma::mat VX = sol.cols(0,q-1);
  arma::mat Vy = sol.col(q);
  arma::mat beta = solve(X.t()*VX,X.t()*Vy);
  arma::mat u = ZK.t()*(Vy-VX*beta);
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*pi/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title Solve Multivariate Model
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Zu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K,
                    double tol=1e-6, int maxIter=1000){
  arma::uword n = Y.n_rows;
  arma::uword m = Y.n_cols;
  arma::mat ZK = Z*K;
  arma::mat ZKZ = ZK*Z.t();
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, ZKZ);
  arma::mat Yt = Y.t()*eigvec;
  arma::mat Xt = X.t()*eigvec;
  arma::mat Vu = cov(Y)/2;
  arma::mat Ve = Vu;
  arma::mat W = Xt.t()*inv_sympd(Xt*Xt.t());
  arma::mat B = Yt*W; //BLUEs
  arma::cube Hinv(m,m,n);
  arma::mat tolEye = tol*arma::eye(m,m);
  arma::mat Gt(m,n), sigma(m,m), BNew,
  VeNew(m,m), VuNew(m,m);
  double denom, numer;
  bool converging=true;
  int iter=0;
  while(converging){
    ++iter;
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    // The same m by m inverse is needed again in the second loop and
    // for the BLUPs, so it is formed once per iteration and kept
    for(arma::uword i=0; i<n; ++i){
      Hinv.slice(i) = inv_sympd(eigval(i)*Vu+Ve+tolEye);
      Gt.col(i) = (eigval(i)*Vu)*Hinv.slice(i)*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(arma::uword i=0; i<n; ++i){
      arma::mat lVu = eigval(i)*Vu;
      sigma = lVu-lVu*Hinv.slice(i)*lVu;
      arma::mat resid = Yt.col(i)-BNew*Xt.col(i)-Gt.col(i);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*(resid*resid.t()+sigma);
    }
    denom = fabs(sum(Ve.diag()));
    if(denom>0.0){
      numer = fabs(sum(VeNew.diag()-Ve.diag()));
      if((numer/denom)<tol) converging=false;
    }
    Ve = VeNew;
    Vu = VuNew;
    B = BNew;
    if(iter>=maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  // The BLUPs follow from the eigendecomposition already computed. The
  // variance of vec(E) is (eigvec (x) I)*blockdiag(eigval_i*Vu+Ve)*
  // (eigvec' (x) I), so its inverse is applied by rotating, solving each
  // m by m block, and rotating back. The n*m by n*m matrix that the
  // Kronecker form would build is never needed.
  for(arma::uword i=0; i<n; ++i){
    Hinv.slice(i) = inv_sympd(eigval(i)*Vu+Ve+tolEye);
  }
  arma::mat E = Y.t() - B*X.t();
  arma::mat Et = E*eigvec;
  arma::mat Wt(m,n);
  for(arma::uword i=0; i<n; ++i){
    Wt.col(i) = Hinv.slice(i)*Et.col(i);
  }
  arma::mat U = Vu*(Wt*eigvec.t())*ZK; //BLUPs
  //Log Likelihood calculation. The quadratic form is a sum over the
  //rotated blocks, and the determinant of the Kronecker sum is the
  //product of the determinants of those blocks.
  double ll = -0.5*accu(Et%Wt) - double(n*m)/2.0*log(2*pi);
  double value;
  double sign;
  for(arma::uword i=0; i<n; ++i){
    log_det(value, sign, eigval(i)*Vu+Ve);
    ll -= 0.5*value*sign;
  }
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("LL")=ll,
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multikernel Model
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Zlist a list of Z matrices
//' @param Klist a list of K matrices
//' @param maxIter maximum number of iteration
//' @param tol tolerance for convergence
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMKM(arma::mat& y, arma::mat& X,
                    arma::field<arma::mat>& Zlist,
                    arma::field<arma::mat>& Klist,
                    int maxIter=40, double tol=1e-4){
  arma::uword k = Klist.n_elem;
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  double df = double(n)-double(q);
  arma::field<arma::mat> V(k);
  for(arma::uword i=0; i<k; ++i){
    V(i) = Zlist(i)*Klist(i)*Zlist(i).t();
  }
  // Choose once between forming the n by n products WQX*V(i) and working
  // from Zlist and Klist. Both give the same average information matrix.
  // The second carries extra products by Klist, so it only pays off when
  // the number of columns is well below n, roughly below 0.4*n.
  double denseCost=0.0, thinCost=0.0;
  for(arma::uword i=0; i<k; ++i){
    double mi = double(Klist(i).n_cols);
    denseCost += double(n)*double(n)*double(n);
    thinCost += double(n)*double(n)*mi + double(n)*mi*mi;
    for(arma::uword j=i; j<k; ++j){
      double mj = double(Klist(j).n_cols);
      denseCost += double(n)*double(n);
      thinCost += double(n)*mi*mj + mi*mi*mj + mi*mj*mj;
    }
  }
  bool useThin = thinCost<denseCost;
  arma::field<arma::mat> C(k);
  arma::mat A(k+1,k+1), W0(n,n), W(n,n), WX(n,q), WQX(n,n);
  arma::vec qvec(k+1), sigma(k+1);
  arma::mat g(n,1);
  double logdetV=0;
  double rss, ldet, llik, llik0=0, deltaLlik, taper,
  value, sign;
  bool invPass;
  arma::field<arma::mat> T(k);
  sigma.fill(var(y.col(0)));
  int iter = 0;
  while(true){
    ++iter;
    W0 = V(0)*sigma(0);
    W0.diag() += sigma(k);
    for(arma::uword i=1; i<k; ++i){
      W0 += V(i)*sigma(i);
    }
    invPass = invAndLogDet(W,logdetV,W0);
    if(!invPass){
      W = pinv(W0);
      log_det(value, sign, W0);
      logdetV = value*sign;
    }
    WX = W*X;
    WQX = W - WX*solve(X.t()*WX, WX.t());
    g = WQX*y;
    rss = as_scalar(y.t()*g);
    sigma = sigma*(rss/df);
    WQX = WQX*(df/rss);
    g = g*(df/rss);
    // WQX has rank n-q, so its ordinary determinant is zero and
    // log_det cannot be used on it. What the likelihood needs is the
    // product of its nonzero eigenvalues, which follows from
    // |WQX|+ = |X'X|/(|V|*|X'V^-1*X|). The constant |X'X| is dropped
    // because only changes in llik are used.
    log_det(value, sign, X.t()*WX);
    ldet = -logdetV - value*sign + df*log(df/rss);
    llik = ldet/2 - df/2;
    if(iter == 1) llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    // y'*T_i*WQX*y equals (WQX*y)'*V_i*(WQX*y), which is quadratic
    // rather than cubic in n. A is symmetric, so only its upper
    // triangle is computed, and WQX is symmetric, so the transposes
    // that would each need a temporary n by n matrix are dropped.
    if(useThin){
      // V(i) = Zlist(i)*Klist(i)*Zlist(i)', so every trace the average
      // information matrix needs is reachable through C(i) = WQX*Zlist(i)
      // and D = Zlist(i).t()*C(j), without forming WQX*V(i):
      //   tr(WQX*V_i)         = accu(Klist(i) % D_ii)
      //   tr(WQX*V_i*WQX*V_j) = accu((D'*Klist(i)*D) % Klist(j))
      //   tr(WQX*V_i*WQX)     = accu(Klist(i) % (C(i).t()*C(i)))
      for(arma::uword i=0; i<k; ++i){
        C(i) = WQX*Zlist(i);
      }
      for(arma::uword i=0; i<k; ++i){
        arma::mat zg = Zlist(i).t()*g;
        arma::mat Dii = Zlist(i).t()*C(i);
        qvec(i) = as_scalar(zg.t()*Klist(i)*zg) - accu(Klist(i)%Dii);
        A(i,i) = accu((Dii.t()*(Klist(i)*Dii))%Klist(i));
        for(arma::uword j=i+1; j<k; ++j){
          arma::mat Dij = Zlist(i).t()*C(j);
          A(i,j) = accu((Dij.t()*(Klist(i)*Dij))%Klist(j));
          A(j,i) = A(i,j);
        }
        A(i,k) = accu(Klist(i)%(C(i).t()*C(i)));
        A(k,i) = A(i,k);
      }
    }else{
      for(arma::uword i=0; i<k; ++i){
        T(i) = WQX*V(i);
      }
      for(arma::uword i=0; i<k; ++i){
        qvec(i) = as_scalar(g.t()*V(i)*g) - sum(T(i).diag());
        for(arma::uword j=i; j<k; ++j){
          A(i,j) = accu(T(i)%T(j).t());
          A(j,i) = A(i,j);
        }
        A(i,k) = accu(T(i)%WQX);
        A(k,i) = A(i,k);
      }
    }
    A(k,k) = accu(WQX%WQX);
    qvec(k) = as_scalar(g.t()*g) - sum(WQX.diag());
    A = pinv(A);
    qvec = A*qvec;
    if(iter == 1){
      taper = 0.5;
    }else if(iter == 2){
      taper = 0.7;
    }else{
      taper = 0.9;
    }
    sigma += taper*qvec;
    while(sigma.min() < -(1e-6)){
      sigma(sigma.index_min()) = -(1e-6);
    }
    if((iter > 1) & (fabs(deltaLlik) < tol*10)){
      break;
    }
    if(max(abs(qvec)) < tol){
      break;
    }
    if(iter >= maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  while(sigma.min() < 0.0){
    sigma(sigma.index_min()) = 0.0;
  }
  arma::mat beta(q,1), ee(n,1);
  arma::field<arma::mat> u(k);
  beta = solve(X.t()*W*X,X.t()*W*y);
  ee = y - X*beta;
  for(arma::uword i=0; i<k; ++i){
    u(i) = (Klist(i)*sigma(i))*Zlist(i).t()*W*ee;
  }
  arma::vec Vu(k), Ve(1);
  Vu = sigma(arma::span(0,k-1));
  Ve = sigma(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=llik,
                            Rcpp::Named("iter")=iter);
}
