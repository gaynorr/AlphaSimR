// solveUVM and solveMVM are based on R/EMMREML functions
// solveMKM is based on the mmer function in R/sommer
#include "alphasimr.h"

#ifdef ARMA_USE_LAPACK

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
#endif


// Replacement for Armadillo's eig_sym
// Fixes an error with decompisition of large matrices
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

// Generates weighted matrix
// Allows for heterogenous variance due to unequal replication
void sweepReps(arma::mat& X, arma::vec reps){
  reps = sqrt(reps);
  for(arma::uword i=0; i<X.n_cols; ++i){
    X.col(i) = X.col(i)%reps;
  }
}

Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M){
  arma::uword n = y.n_rows;
  arma::uword q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));

  // Construct system of equations for eigendecomposition
  arma::mat S = -(X*inv_sympd(X.t()*X)*X.t());
  S.diag() += 1;
  arma::mat H = M*M.t(); // Used later
  H.diag() += offset;
  S = S*H*S;

  // Compute eigendecomposition
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, S);

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
  H = inv_sympd(H);
  arma::mat XH = X.t()*H;
  arma::mat beta = solve(XH*X,XH*y);
  arma::mat u = M.t()*(H*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u);
}

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
  arma::mat Gt(m,n), sigma(m,m), BNew, 
  VeNew(m,m), VuNew(m,m);
  double denom, numer;
  bool converging=true;
  int iter=0;
  while(converging){
    ++iter;
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    for(arma::uword i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(arma::uword i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(eigval(i)*Vu);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*((Yt.col(i)-BNew*Xt.col(i)-Gt.col(i))*
        (Yt.col(i)-BNew*Xt.col(i)-Gt.col(i)).t()+sigma);
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
  arma::mat HI = inv_sympd(kron(M*M.t(), Vu)+
    kron(arma::eye(n,n), Ve)+
    tol*arma::eye(n*m,n*m));
  arma::mat E = Y.t() - B*X.t();
  arma::mat U = kron(arma::eye(M.n_cols,M.n_cols), Vu)*kron(M.t(),
                     arma::eye(m,m))*(HI*vectorise(E)); //BLUPs
  U.reshape(m,U.n_elem/m);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("iter")=iter);
}

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
  arma::mat A(k+1,k+1), W0(n,n), W(n,n), WX(n,q), WQX(n,n);
  arma::vec qvec(k+1), sigma(k+1);
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
    invPass = inv_sympd(W,W0);
    if(!invPass){
      W = pinv(W0);
    }
    WX = W*X;
    WQX = W - WX*solve(X.t()*WX, WX.t());
    rss = as_scalar(y.t()*WQX*y);
    sigma = sigma*(rss/df);
    WQX = WQX*(df/rss);
    log_det(value, sign, WQX);
    ldet = value*sign;
    llik = ldet/2 - df/2;
    if(iter == 1) llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    for(arma::uword i=0; i<k; ++i){
      T(i) = WQX*V(i);
    }
    for(arma::uword i=0; i<k; ++i){
      qvec(i) = as_scalar(y.t()*T(i)*WQX*y - sum(T(i).diag()));
      for(arma::uword j=0; j<k; ++j){
        A(i,j) = accu(T(i)%T(j).t());
      }
      A(i,k) = accu(T(i)%WQX.t());
    }
    for(arma::uword j=0; j<k; ++j){
      A(k,j) = accu(WQX%T(j).t());
    }
    A(k,k) = accu(WQX%WQX.t());
    qvec(k) = as_scalar(y.t()*WQX*WQX*y - sum(WQX.diag()));
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

//Uses EM algorithm to solve a mixed model with 1 random effect
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
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m-1)) = X.t()*M;
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) = M.t()*X;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M.t()*M;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M.t()*Y;
  arma::mat RHSinv = inv(RHS);
  LHS = RHSinv*Rvec;
  if(useEM){
    VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
    VuN = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
    )/m;
    delta = VeN/VuN-lambda;
    while(fabs(delta)>tol){
      Ve = VeN;
      Vu = VuN;
      RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += delta;
      lambda += delta;
      RHSinv = inv(RHS);
      LHS = RHSinv*Rvec;
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
      VuN = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
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

//Uses EM algorithm to solve a mixed model with 2 random effects
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
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) = M1.t()*X;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M1.t()*M1;
  RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)) = M1.t()*M2;
  // Third row
  RHS(arma::span(q+m,q+2*m-1),arma::span(0,q-1)) = M2.t()*X;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q,q+m-1)) = M2.t()*M1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)) = M2.t()*M2;
  // Add to diagonal
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += lambda2;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M1.t()*Y;
  Rvec(arma::span(q+m,q+2*m-1),0) = M2.t()*Y;
  arma::mat RHSinv = inv(RHS);
  LHS = RHSinv*Rvec;
  if(useEM){
    VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
    Vu1N = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
    )/m;
    Vu2N = as_scalar(
      LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
        Ve*sum(RHSinv(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag())
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
      RHSinv = inv(RHS);
      LHS = RHSinv*Rvec;
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
      Vu1N = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
      )/m;
      Vu2N = as_scalar(
        LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
          Ve*sum(RHSinv(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag())
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

//Uses EM algorithm to solve a mixed model with 3 random effects
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
    arma::mat u(m,2);
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
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) = M1.t()*X;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M1.t()*M1;
  RHS(arma::span(q,q+m-1),arma::span(q+m,q+2*m-1)) = M1.t()*M2;
  RHS(arma::span(q,q+m-1),arma::span(q+2*m,q+3*m-1)) = M1.t()*M3;
  // Third row
  RHS(arma::span(q+m,q+2*m-1),arma::span(0,q-1)) = M2.t()*X;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q,q+m-1)) = M2.t()*M1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)) = M2.t()*M2;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+2*m,q+3*m-1)) = M2.t()*M3;
  // Fourth row
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(0,q-1)) = M3.t()*X;
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q,q+m-1)) = M3.t()*M1;
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+m,q+2*m-1)) = M3.t()*M2;
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)) = M3.t()*M3;
  // Add to diagonal
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda1;
  RHS(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag() += lambda2;
  RHS(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)).diag() += lambda3;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M1.t()*Y;
  Rvec(arma::span(q+m,q+2*m-1),0) = M2.t()*Y;
  Rvec(arma::span(q+2*m,q+3*m-1),0) = M3.t()*Y;
  arma::mat RHSinv = inv(RHS);
  LHS = RHSinv*Rvec;
  if(useEM){
    VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
    Vu1N = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
    )/m;
    Vu2N = as_scalar(
      LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
        Ve*sum(RHSinv(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag())
    )/m;
    Vu3N = as_scalar(
      LHS(arma::span(q+2*m,q+3*m-1),0).t()*LHS(arma::span(q+2*m,q+3*m-1),0)+
        Ve*sum(RHSinv(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)).diag())
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
      RHSinv = inv(RHS);
      LHS = RHSinv*Rvec;
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
      Vu1N = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
      )/m;
      Vu2N = as_scalar(
        LHS(arma::span(q+m,q+2*m-1),0).t()*LHS(arma::span(q+m,q+2*m-1),0)+
          Ve*sum(RHSinv(arma::span(q+m,q+2*m-1),arma::span(q+m,q+2*m-1)).diag())
      )/m;
      Vu3N = as_scalar(
        LHS(arma::span(q+2*m,q+3*m-1),0).t()*LHS(arma::span(q+2*m,q+3*m-1),0)+
          Ve*sum(RHSinv(arma::span(q+2*m,q+3*m-1),arma::span(q+2*m,q+3*m-1)).diag())
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

// Called by fastRRBLUP function
// An implementation of the Gauss-Seidel method for solving 
// mixed model equations for an RR-BLUP model
// [[Rcpp::export]]
Rcpp::List callFastRRBLUP(arma::vec y,
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                          double Vu, double Ve, arma::uword maxIter, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::mat Md(y.n_rows,1);
  arma::rowvec Mdr(M.n_cols);
  arma::rowvec Mmean(M.n_cols);
  arma::vec X(y.n_rows);
  double lhs, rhs, eps, beta=0, solOld;
  arma::rowvec XpX(M.n_cols);
  for(arma::uword i=0; i<M.n_cols; ++i){
    Md = genoToGenoA(M.col(i), ploidy, 1);
    Mmean(i) = as_scalar(mean(Md));
    Md -= Mmean(i);
    XpX(i) = accu(Md%Md)/Ve;
  }
  double lambda = 1/Vu;
  arma::vec u(M.n_cols);
  u.fill(1e-6);
  arma::vec e = y;
  double OpO = M.n_rows/Ve;
  arma::uvec order = arma::regspace<arma::uvec>(0,M.n_cols-1);
  arma::uword k;
  for(arma::uword iter=0; iter<maxIter; ++iter){
    e += beta;
    rhs = accu(e)/Ve;
    beta = rhs/OpO;
    e -= beta;
    eps=0;
    order = shuffle(order);
    for(arma::uword i=0; i<M.n_cols; ++i){
      k = order(i);
      Md = genoToGenoA(M.col(k), ploidy, 1);
      Md -= Mmean(k);
      e += Md*u(k);
      lhs = XpX(k)+lambda;
      rhs = accu(Md%e/Ve);
      solOld = u(k);
      u(k) = rhs/lhs;
      e -= Md*u(k);
      eps += pow((u(k)-solOld),2);
    }
    if(iter%200 ==0){
      for(arma::uword i=0; i<M.n_rows; ++i){
        Mdr = genoToGenoA(M.row(i), ploidy, nThreads);
        Mdr -= Mmean;
        X(i) = as_scalar(Mdr*u);
      }
      e = y-X-beta;
    }
    if(eps<1e-8){
      break;
    }
  }
  return Rcpp::List::create(Rcpp::Named("alpha")=u,
                            Rcpp::Named("beta")=-as_scalar(Mmean*u),
                            Rcpp::Named("mu")=beta);
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x, arma::vec reps,
                      arma::field<arma::Cube<unsigned char> >& geno, 
                      arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                      bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  if(useReps){
    sweepReps(y,reps);
    sweepReps(X,reps);
    sweepReps(M,reps);
  }
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
Rcpp::List callRRBLUP2(arma::mat y, arma::uvec x, arma::vec reps,
                       arma::field<arma::Cube<unsigned char> >& geno, 
                       arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                       double Vu, double Ve, double tol, int maxIter, 
                       bool useEM, bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  if(useReps){
    sweepReps(y,reps);
    sweepReps(X,reps);
    sweepReps(M,reps);
  }
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
Rcpp::List callRRBLUP_D(arma::mat y, arma::uvec x, arma::vec reps,
                        arma::field<arma::Cube<unsigned char> >& geno, 
                        arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                        int maxIter, bool useReps, int nThreads){
  // Fit GS model
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Mlist(0));
  Mlist(1) = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Mlist(1),1));
  if(useReps){
    sweepReps(y,reps);
    sweepReps(X,reps);
    sweepReps(Mlist(0),reps);
    sweepReps(Mlist(1),reps);
  }
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
Rcpp::List callRRBLUP_D2(arma::mat y, arma::uvec x, arma::vec reps,
                         arma::field<arma::Cube<unsigned char> >& geno, 
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc,
                         int maxIter, double Va, double Vd, double Ve, 
                         double tol, bool useEM, bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::Mat<unsigned char> M = getGeno(geno,lociPerChr,lociLoc,nThreads);
  arma::mat Ma = genoToGenoA(M,ploidy,nThreads);
  arma::rowvec Mmean = mean(Ma);
  arma::mat Md = genoToGenoD(M,ploidy,nThreads);
  arma::mat X;
  X = join_rows(makeX(x),sum(Md,1));
  if(useReps){
    sweepReps(y,reps);
    sweepReps(X,reps);
    sweepReps(Ma,reps);
    sweepReps(Md,reps);
  }
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
Rcpp::List callRRBLUP_MV(arma::mat Y, arma::uvec x, arma::vec reps,
                         arma::field<arma::Cube<unsigned char> >& geno, 
                         arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                         int maxIter, bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat X = makeX(x);
  arma::mat M = genoToGenoA(getGeno(geno,lociPerChr,lociLoc,nThreads),
                            ploidy,nThreads);
  arma::rowvec Mmean = mean(M);
  if(useReps){
    sweepReps(Y,reps);
    sweepReps(X,reps);
    sweepReps(M,reps);
  }
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
Rcpp::List callRRBLUP_GCA(arma::mat y, arma::uvec x, arma::vec reps,
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                          int maxIter, bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::field<arma::mat> Mlist(2);
  Mlist(0) = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  Mlist(1) = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
                         nThreads),ploidy/2,nThreads);
  arma::rowvec Mmean1 = mean(Mlist(0));
  arma::rowvec Mmean2 = mean(Mlist(1));
  arma::mat X = makeX(x);
  if(useReps){
    sweepReps(y, reps);
    sweepReps(X, reps);
    sweepReps(Mlist(0), reps);
    sweepReps(Mlist(1), reps);
  }
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
Rcpp::List callRRBLUP_GCA2(arma::mat y, arma::uvec x, arma::vec reps,
                           arma::field<arma::Cube<unsigned char> >& geno, 
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                           int maxIter, double Vu1, double Vu2, double Ve, 
                           double tol, bool useEM, bool useReps, int nThreads){
  arma::uword ploidy = geno(0).n_cols;
  arma::mat M1 = genoToGenoA(getMaternalGeno(geno,lociPerChr,lociLoc,
                                             nThreads),ploidy/2,nThreads);
  arma::mat M2 = genoToGenoA(getPaternalGeno(geno,lociPerChr,lociLoc,
                                             nThreads),ploidy/2,nThreads);
  arma::rowvec Mmean1 = mean(M1);
  arma::rowvec Mmean2 = mean(M2);
  arma::mat X = makeX(x);
  if(useReps){
    sweepReps(y, reps);
    sweepReps(X, reps);
    sweepReps(M1, reps);
    sweepReps(M2, reps);
  }
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
Rcpp::List callRRBLUP_SCA(arma::mat y, arma::uvec x, arma::vec reps,
                          arma::field<arma::Cube<unsigned char> >& geno, 
                          arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                          int maxIter, bool useReps, int nThreads){
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
  if(useReps){
    sweepReps(y, reps);
    sweepReps(X, reps);
    sweepReps(Mlist(0), reps);
    sweepReps(Mlist(1), reps);
    sweepReps(Mlist(2), reps);
  }
  
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
Rcpp::List callRRBLUP_SCA2(arma::mat y, arma::uvec x, arma::vec reps,
                           arma::field<arma::Cube<unsigned char> >& geno, 
                           arma::Col<int>& lociPerChr, arma::uvec lociLoc, 
                           int maxIter, double Vu1, double Vu2, double Vu3, 
                           double Ve, double tol, bool useEM, bool useReps, 
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
  if(useReps){
    sweepReps(y, reps);
    sweepReps(X, reps);
    sweepReps(M1, reps);
    sweepReps(M2, reps);
    sweepReps(M3, reps);
  }
  
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
