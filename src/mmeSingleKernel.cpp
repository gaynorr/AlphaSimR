#include "mme.h"

// Single trait, single kernel mixed model with a known covariance among
// levels. It is written from the published descriptions of FaST-LMM
// (Lippert et al. 2011, doi:10.1038/nmeth.1681) and of restricted maximum
// likelihood (Patterson and Thompson 1971, doi:10.1093/biomet/58.3.545), and
// shares its likelihood machinery with solveRRBLUP, in MME.cpp, through
// mme.h.

//' @title Solve Univariate Model
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Zu+e}, where
//' \eqn{u \sim N(0,K\sigma^2_u)} and \eqn{e \sim N(0,I\sigma^2_e)}. The
//' variance components are estimated by REML
//' \insertCite{patterson1971Recovery}{AlphaSimR} using the spectral approach
//' of FaST-LMM \insertCite{lippert2011FaST}{AlphaSimR}, applied to
//' \eqn{ZKZ'}. The fixed effects are the generalized least squares estimates
//' and the random effects are the BLUPs \insertCite{henderson1975Best}{AlphaSimR},
//' both at the estimated variance components.
//'
//' @details
//' The decomposition is taken in whichever space is smaller. When there are
//' fewer levels than records, \eqn{ZKZ'} is never formed: its nonzero
//' eigenvalues are found from \eqn{K} and \eqn{Z'Z}, and the directions in
//' which it is zero are summarized analytically. \eqn{K} may be singular,
//' and records may repeat levels.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns. Columns that are linearly
//' dependent on the others are dropped and get a fixed effect of zero. An
//' intercept, a column whose values are all the same, is always kept.
//' @param Z a matrix with n rows and m columns
//' @param K a positive semidefinite matrix with m rows and m columns.
//' Eigenvalues that are negative only at rounding level are treated as zero.
//'
//' @return A list with the random effect variance \code{Vu}, the residual
//' variance \code{Ve}, the fixed effects \code{beta}, the BLUPs \code{u} of
//' the m levels, and the restricted log-likelihood \code{LL}. \code{LL} is
//' \deqn{-\frac{1}{2}\left[(n-p)\log(2\pi)+\log|V|+\log|X'V^{-1}X|+y'Py\right]}
//' for \eqn{V=ZKZ'\sigma^2_u+I\sigma^2_e}, \eqn{p} the rank of \eqn{X} and
//' \eqn{P=V^{-1}-V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}}, with \eqn{X} reduced to
//' the columns that are kept. It excludes the \eqn{\log|X'X|} term that some
//' definitions include.
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K){
  arma::uword n = y.n_rows;
  arma::uword m = Z.n_cols;
  if( (X.n_rows!=n) || (Z.n_rows!=n) || (K.n_rows!=m) || (K.n_cols!=m) ){
    Rcpp::stop("solveUVM: the dimensions of y, X, Z and K do not agree");
  }
  // Duplicated or collinear columns leave X'*H^-1*X singular, or singular
  // only up to rounding, in which case its Cholesky factorization can succeed
  // and return meaningless fixed effects. Such columns are therefore dropped
  // before fitting, and everything below works with Xr, which spans the same
  // space with full column rank.
  const arma::uvec xKeep = independentColumns(X);
  if(xKeep.n_elem==0){
    Rcpp::stop("solveUVM: X must have at least one nonzero column");
  }
  const arma::mat Xr = X.cols(xKeep);
  arma::uword q = Xr.n_cols;
  if(q>=n){
    Rcpp::stop("solveUVM: there must be more records than fixed effects");
  }
  double df = double(n)-double(q);

  arma::vec yv = y.col(0);

  arma::vec s2;  // eigenvalues of Z*K*Z' that are kept
  arma::vec yt;  // U'*y
  arma::mat Xt;  // U'*X
  arma::mat U;   // eigenvectors of Z*K*Z', full rank branch only
  arma::mat W;   // Z*K^(1/2), low rank branch only
  arma::mat WtU; // W'*U, low rank branch only
  arma::mat KhalfQ; // K^(1/2) as Q*diag(sqrt(lambda)), low rank branch only
  // Quantities outside the span of the kept eigenvectors
  double yy0 = 0.0;
  arma::vec Xy0(q, arma::fill::zeros);
  arma::mat XX0(q, q, arma::fill::zeros);
  double nNull = 0.0;

  if(m<n){
    // Writing K = Q*diag(lambda)*Q' gives Z*K*Z' = W*W' for
    // W = Z*Q*diag(sqrt(lambda)), and the model becomes an RR-BLUP model in
    // W. Only the k positive eigenvalues of K contribute, so W has k columns
    // and the rank of Z*K*Z' is at most k. Its nonzero eigenvalues are those
    // of W'*W, which is k by k, and the eigenvectors in R^n follow as
    // W*V*diag(1/s), so nothing n by n is ever formed.
    arma::vec lambda(m);
    arma::mat Q(m,m);
    if(eigen2(lambda, Q, K)!=0){
      Rcpp::stop("solveUVM: eigendecomposition of K failed");
    }
    double maxLambda = arma::abs(lambda).max();
    if(lambda.min() < -sqrt(arma::datum::eps)*maxLambda){
      Rcpp::stop("solveUVM: K is not positive semidefinite");
    }
    arma::uvec keepK = find(lambda>(double(m)*arma::datum::eps*maxLambda));
    if(keepK.n_elem==0){
      Rcpp::stop("solveUVM: K has no positive eigenvalues");
    }
    KhalfQ = Q.cols(keepK);
    KhalfQ.each_row() %= sqrt(lambda(keepK)).t();
    Q.reset();
    W = Z*KhalfQ;
    arma::uword k = W.n_cols;

    arma::vec Wty = W.t()*yv;
    arma::mat WtX = W.t()*Xr;
    arma::vec eigval(k);
    arma::mat eigvec(k,k);
    if(eigen2(eigval, eigvec, W.t()*W)!=0){
      Rcpp::stop("solveUVM: eigendecomposition of Z*K*Z' failed");
    }
    double maxEig = eigval.max();
    if(maxEig<0.0) maxEig = 0.0;
    arma::uvec keep = find(eigval>(double(k)*arma::datum::eps*maxEig));
    if(keep.n_elem==0){
      Rcpp::stop("solveUVM: Z*K*Z' is zero");
    }
    s2 = eigval(keep);
    arma::mat Vk = eigvec.cols(keep);
    arma::vec sInv = 1.0/sqrt(s2);
    yt = sInv%(Vk.t()*Wty);
    Xt = Vk.t()*WtX;
    Xt.each_col() %= sInv;
    WtU = Vk;
    WtU.each_row() %= sqrt(s2).t();
    // The remaining n-k directions have a zero eigenvalue. Their
    // contribution is the same for every delta apart from a factor of
    // 1/delta, so it is summarized here rather than being stored.
    nNull = double(n)-double(s2.n_elem);
    yy0 = accu(yv%yv)-accu(yt%yt);
    if(yy0<0.0) yy0 = 0.0;
    Xy0 = Xr.t()*yv-Xt.t()*yt;
    XX0 = Xr.t()*Xr-Xt.t()*Xt;
    XX0 = 0.5*(XX0+XX0.t()); // rounding only
  }else{
    // Z*K*Z' is no larger than K, so it is decomposed directly. Its
    // eigenvectors span all of R^n, so there is no null space to account
    // for, and zero eigenvalues are harmless because delta is positive.
    arma::mat G = kernelOnRecords(Z, K);
    G = 0.5*(G+G.t()); // rounding only
    arma::vec eigval(n);
    U.set_size(n,n);
    if(eigen2(eigval, U, G)!=0){
      Rcpp::stop("solveUVM: eigendecomposition of Z*K*Z' failed");
    }
    G.reset();
    double maxEig = arma::abs(eigval).max();
    if(eigval.min() < -sqrt(arma::datum::eps)*maxEig){
      Rcpp::stop("solveUVM: K is not positive semidefinite");
    }
    eigval.elem(find(eigval<0.0)).zeros(); // rounding only
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

  // The residual variance is profiled out analytically, so REML reduces to
  // a one dimensional search over delta = Ve/Vu
  double delta = fastLMMDelta(args);

  FastLMMFit fit = fastLMMFit(delta, yt, Xt, s2, yy0, Xy0, XX0, nNull);
  if(!fit.ok){
    Rcpp::stop("solveUVM: mixed model equations could not be solved");
  }
  arma::mat beta = expandFixedEffects(fit.beta, xKeep, X.n_cols);

  // The BLUPs are u = K*Z'*H^-1*(y-X*beta), where V = Vu*H, which the
  // decomposition already in hand applies without a second factorization
  arma::vec et = yt-Xt*fit.beta; // U'*(y-X*beta)
  arma::vec etd = et/(s2+delta);
  arma::mat u;
  if(m<n){
    // K*Z' = KhalfQ*W', so u = KhalfQ*a for the RR-BLUP effects a = W'*H^-1*e
    arma::vec Wte = W.t()*(yv-Xr*fit.beta);
    arma::vec a = WtU*etd+(Wte-WtU*et)/delta;
    u = KhalfQ*a;
  }else{
    u = K*(Z.t()*(U*etd));
  }

  // At the optimum Vu = r/(n-p), which makes y'Py equal to n-p and leaves
  // the restricted log-likelihood in closed form
  double Vu = fit.r/df;
  double Ve = delta*Vu;
  double LL = -0.5*(df*(log(2.0*arma::datum::pi)+log(Vu)+1.0) +
                    fit.logDetH + fit.logDetA);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=LL);
}
