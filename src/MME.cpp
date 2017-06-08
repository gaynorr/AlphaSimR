// solveUVM and solveMVM are based on R/EMMREML functions
// solveMKM is based on the mmer function in R/sommer 
#include "alphasimr.h"
#include <iostream>
#include <string>

// Note: Fortran compiler appends '_' to subroutine name
// See http://www.netlib.org/lapack/explore-html/ for description of args
extern "C" void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL,
                       double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z,
                       int* LDZ, int* ISUPPZ, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

// Replacement for Armadillo's eig_sym
// Fixes an errors with decompisition of large matrices on Eddie
// If calcVec = false, eigvec is not used
// It would be better to template this function
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
  int N = X.n_rows;
  // A = X
  int LDA = N;
  double VL = 0.0;
  double VU = 0.0;
  int IL;
  int IU;
  double ABSTOL = 0.0;
  int M = N;
  // W=eigval
  // Z=eigvec
  int LDZ = N;
  arma::ivec ISUPPZ(2*M);
  // WORK length to be determined
  double tmpWORK;
  int LWORK = -1; // To be calculated
  // IWORK length to be determined
  int tmpIWORK;
  int LIWORK = -1; // To be calculated
  int INFO;
  // Calculate LWORK and LIWORK
  dsyevr_(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
          &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&tmpWORK,&LWORK,&tmpIWORK,&LIWORK,&INFO);
  LWORK = int(tmpWORK);
  LIWORK = tmpIWORK;
  // Allocate WORK and IWORK
  arma::vec WORK(LWORK);
  arma::ivec IWORK(LIWORK);
  // Perform decomposition
  dsyevr_(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
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

/*
 * Reads a text file into an arma::mat
 * Must be supplied with the correct number of rows and columns for the matrix
 * A header can be skipped by setting skipRows to 1
 * Row names can be skipped by setting skipCols to 1
 */
arma::mat readMat(std::string fileName, int rows, int cols, 
                  char sep=' ', int skipRows=0, int skipCols=0){
  arma::mat output(rows,cols);
  std::ifstream file(fileName.c_str());
  std::string line;
  //Skip rows
  for(int i=0; i<skipRows; ++i){
    std::getline(file,line);
  }
  //Read rows
  for(int i=0; i<rows; ++i){
    std::getline(file,line);
    std::stringstream lineStream(line);
    std::string cell;
    //Skip columns
    for(int j=0; j<skipCols; ++j){
      std::getline(lineStream,cell,sep);
    }
    //Read columns
    for(int j=0; j<cols; ++j){
      std::getline(lineStream,cell,sep);
      output(i,j) = std::atof(cell.c_str());
      //output(i,j) = std::stod(cell);
    }
  }
  file.close();
  return output;
}

// Produces a sum to zero design matrix with an intercept
arma::mat makeX(arma::uvec& x){
  int nTrain = x.n_elem;
  double nLevels = x.max();
  arma::mat X(nTrain,nLevels);
  if(nLevels==1){
    X.ones();
  }else{
    X.zeros();
    X.col(0).ones();
    for(int i=0; i<nTrain; ++i){
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
arma::mat makeZ(arma::uvec& z, int nGeno){
  int nTrain = z.n_elem;
  arma::mat Z(nTrain,nGeno,arma::fill::zeros);
  for(int i=0; i<nTrain; ++i){
    Z(i,z(i)) = 1;
  }
  return Z;
}

// Generates weighted matrix
// Allows for heterogenous variance due to unequal replication
void sweepReps(arma::mat& X, arma::vec& reps){
  for(int i=0; i<X.n_cols; ++i){
    X.col(i) = X.col(i)/reps;
  }
}

//' @title Solve Univariate Model
//' 
//' @description
//' Solves a univariate mixed model of form \deqn{y=X\beta+Zu+e}.
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
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));
  
  // Construct system of equations for eigendecomposition
  arma::mat S = arma::eye(n,n) - X*arma::inv_sympd(X.t()*X)*X.t();
  arma::mat ZK = Z*K;
  arma::mat ZKZ = ZK*Z.t();
  S = S*(ZKZ+offset*arma::eye(n,n))*S;
  
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
  arma::mat Hinv = arma::inv_sympd(ZKZ+delta*arma::eye(n,n));
  arma::mat XHinv = X.t()*Hinv;
  arma::mat beta = solve(XHinv*X,XHinv*y);
  arma::mat u = ZK.t()*(Hinv*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title Solve Multivariate Model
//' 
//' @description
//' Solves a multivariate mixed model of form \deqn{Y=X\beta+Zu+e}.
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//' @param tol tolerance for convergence
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X, 
                    const arma::mat& Z, const arma::mat& K,
                    double tol=1e-6){
  int n = Y.n_rows;
  int m = Y.n_cols;
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
  arma::mat Gt(m,n);
  arma::mat sigma(m,m);
  arma::mat BNew;
  arma::mat VeNew(m,m);
  arma::mat VuNew(m,m);
  double denom;
  double numer;
  bool converging=true;
  while(converging){
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    for(int i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*arma::inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(int i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*arma::inv_sympd(eigval(i)*Vu+
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
  }
  arma::mat HI = inv_sympd(kron(ZKZ, Vu)+kron(arma::eye(n,n), Ve)+
    tol*arma::eye(n*m,n*m));
  arma::mat E = Y.t() - B*X.t();
  arma::mat U = kron(K, Vu)*kron(Z.t(), 
                     arma::eye(m,m))*(HI*arma::vectorise(E)); //BLUPs
  U.reshape(m,U.n_elem/m);
  //Log Likelihood calculation
  arma::mat ll = -0.5*arma::vectorise(E).t()*HI*arma::vectorise(E);
  ll -= double(n*m)/2.0*log(2*PI);
  double value;
  double sign;
  log_det(value, sign, kron(ZKZ, Vu)+kron(arma::eye(n,n), Ve));
  ll -= 0.5*value*sign;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("LL")=double(ll(0,0)));
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
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMKM(arma::mat& y, arma::mat& X, 
                     arma::field<arma::mat>& Zlist, 
                     arma::field<arma::mat>& Klist){
  int maxcyc = 20;
  double tol = 1e-4;
  int k = Klist.n_elem;
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  arma::field<arma::mat> V(k+1);
  for(int i=0; i<k; ++i){
    V(i) = Zlist(i)*Klist(i)*Zlist(i).t();
  }
  V(k) = arma::eye(n,n);
  k += 1;
  arma::mat A(k,k);
  arma::vec qvec(k);
  arma::vec sigma(k);
  arma::mat W(n,n);
  arma::mat WX(n,q);
  arma::mat WQX(n,n);
  arma::mat rss(1,1);
  arma::vec eigval;
  arma::mat tmpMat(1,1);
  double ldet;
  double llik;
  double llik0;
  double deltaLlik;
  double taper;
  arma::field<arma::mat> T(k);
  sigma.fill(var(y.col(0)));
  for(int cycle=0; cycle<maxcyc; ++cycle){
    W = V(0)*sigma(0);
    for(int i=1; i<k; ++i){
      W += V(i)*sigma(i);
    }
    W = arma::inv_sympd(W);
    WX = W*X;
    WQX = W - WX*arma::solve(X.t()*WX, WX.t());
    rss = y.t()*WQX*y;
    for(int i=0; i<k; ++i){
      sigma(i) = sigma(i)*(rss(0,0)/df);
    }
    WQX = WQX*(df/rss(0,0));
    // Compute eigendecomposition
    eigval.set_size(n);
    eigen2(eigval, tmpMat, WQX, false);
    eigval = eigval(arma::span(q,n-1));
    if(eigval(0)<0.0){
      WQX += (tol-eigval(0))*arma::eye(n,n);
      eigval += (tol-eigval(0));
    }
    ldet = accu(log(eigval));
    llik = ldet/2 - df/2;
    if(cycle == 1) 
      llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    for(int i=0; i<k; ++i){
      T(i) = WQX*V(i);
    }
    for(int i=0; i<k; ++i){
      tmpMat = y.t()*T(i)*WQX*y - sum(T(i).diag());
      qvec(i) = tmpMat(0,0);
      for(int j=0; j<k; ++j){
        A(i,j) = accu(T(i)%T(j).t());
      }
    }
    A = pinv(A);
    qvec = A*qvec;
    if(cycle == 1){
      taper = 0.5;
    }else if(cycle == 2){
      taper = 0.7;
    }else{
      taper = 0.9;
    }
    sigma += taper*qvec;
    while(sigma.min() < -(1e-6)){
      sigma(sigma.index_min()) = -(1e-6);
    }
    if(cycle > 1 & fabs(deltaLlik) < tol*10){
      break;
    }
    if(max(abs(qvec)) < tol){
      break;
    }
  }
  arma::mat beta(q,1);
  arma::field<arma::mat> u(k-1);
  arma::mat ee(n,1);
  beta = solve(X.t()*W*X,X.t()*W*y);
  ee = y - X*beta;
  for(int i=0; i<(k-1); ++i){
    u(i) = (Klist(i)*sigma(i))*Zlist(i).t()*W*ee;
  }
  arma::vec Vu(k-1);
  Vu = sigma(arma::span(0,k-2));
  arma::vec Ve(1);
  Ve = sigma(k-1);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=llik);
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x, arma::vec reps, 
                         std::string genoTrain, int nMarker){
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,y.n_elem,nMarker,' ',0,0);
  reps = sqrt(1.0/reps);
  sweepReps(y,reps);
  sweepReps(X,reps);
  sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  return solveUVM(y, X, Z, K);
}

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_MV(arma::mat Y, arma::uvec x, arma::vec reps, 
                            std::string genoTrain, int nMarker){
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,Y.n_rows,nMarker,' ',0,0);
  reps = sqrt(1.0/reps);
  sweepReps(Y,reps);
  sweepReps(X,reps);
  sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  return solveMVM(Y, X, Z, K);
}

// Called by RRBLUP_GCA function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_GCA(arma::mat y, arma::uvec x, arma::vec reps,
                          std::string genoFemale, std::string genoMale, 
                          int nMarker){
  int n = y.n_rows;
  arma::mat X = makeX(x);
  arma::field<arma::mat> Zlist(2);
  Zlist(0) = readMat(genoFemale,n,nMarker,' ',0,0);
  Zlist(1) = readMat(genoMale,n,nMarker,' ',0,0);
  arma::field<arma::mat> Klist(2);
  Klist(0) = arma::eye(nMarker,nMarker);
  Klist(1) = arma::eye(nMarker,nMarker);
  reps = sqrt(1.0/reps);
  sweepReps(y, reps);
  sweepReps(X, reps);
  sweepReps(Zlist(0), reps);
  sweepReps(Zlist(1), reps);
  return solveMKM(y,X,Zlist,Klist);
}

// Called by RRBLUP_SCA function
// [[Rcpp::export]]
Rcpp::List callRRBLUP_SCA(arma::mat y, arma::uvec x, arma::vec reps,
                          std::string genoFemale, std::string genoMale, 
                          int nMarker){
  int n = y.n_rows;
  arma::mat X = makeX(x);
  arma::field<arma::mat> Zlist(3);
  Zlist(0) = readMat(genoFemale,n,nMarker,' ',0,0);
  Zlist(0) = Zlist(0)*2-1;
  Zlist(1) = readMat(genoMale,n,nMarker,' ',0,0);
  Zlist(1) = Zlist(1)*2-1;
  Zlist(2) = Zlist(0)%Zlist(1);
  arma::field<arma::mat> Klist(3);
  Klist(0) = arma::eye(nMarker,nMarker);
  Klist(1) = arma::eye(nMarker,nMarker);
  Klist(2) = arma::eye(nMarker,nMarker);
  reps = sqrt(1.0/reps);
  sweepReps(y, reps);
  sweepReps(X, reps);
  sweepReps(Zlist(0), reps);
  sweepReps(Zlist(1), reps);
  sweepReps(Zlist(2), reps);
  return solveMKM(y,X,Zlist,Klist);
}

//' @title Calculate G Matrix
//' 
//' @description
//' Calculates the genomic relationship matrix.
//'
//' @param X a matrix of marker genotypes scored as 0,1,2
//' 
//' @return a matrix of the realized genomic relationship.
//'
//' @export
// [[Rcpp::export]]
arma::mat calcG(arma::mat X){
  arma::rowvec p = mean(X,0)/2.0;
  X.each_row() -= 2*p;
  arma::mat G = X*X.t();
  G = G/(2.0*sum(p%(1-p)));
  return G;
}

// Calculates a distance matrix from a marker matrix
// Uses binomial theorem trick
// Inspired by code from:
// http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
// First described here: 
// http://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
//' @title Calculate Euclidean distance
//' 
//' @description
//' Calculates a Euclidean distance matrix using a binomial 
//' theorem trick. Results in much faster computation than the 
//' \code{dist} function in package \code{stats}.
//'
//' @param X a numeric matrix
//' 
//' @return a matrix of columnwise distances
//'
//' @export
// [[Rcpp::export]]
arma::mat fastDist(const arma::mat& X){
  arma::colvec Xn =  sum(square(X),1);
  arma::mat D = -2*(X*X.t());
  D.each_col() += Xn;
  D.each_row() += Xn.t();
  D = sqrt(D);
  D.diag().zeros(); //Removes NaN values
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0); //Assuming there won't be any Inf values
  }
  return D; 
}

//' @title Calculate Paired Euclidean distance
//' 
//' @description
//' Calculates a Euclidean distance between two matrices using 
//' a binomial theorem trick. 
//'
//' @param X a numeric matrix
//' @param Y a numeric matrix
//' 
//' @return a matrix of columnwise distances between matrices 
//' X and Y
//'
//' @export
// [[Rcpp::export]]
arma::mat fastPairDist(const arma::mat& X, const arma::mat& Y){
  arma::colvec Xn =  sum(square(X),1);
  arma::colvec Yn =  sum(square(Y),1);
  arma::mat D = -2*(X*Y.t());
  D.each_col() += Xn;
  D.each_row() += Yn.t();
  D = sqrt(D);
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0); //Assuming there won't be any Inf values
  }
  return D; 
}

//' @title Calculate Gaussian Kernel
//' 
//' @description
//' Calculates a Gaussian kernel using a Euclidean distance 
//' matrix.
//'
//' @param D a matrix of Euclidean distances, 
//' see \code{\link{fastDist}}
//' @param theta the tuning parameter
//' 
//' @return a numeric matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat gaussKernel(arma::mat& D, double theta){
  return exp(-1.0*square(D/theta));
}