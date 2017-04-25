/*
 * Mixed model code from:
 * "A two-part strategy for using genomic selection to develop inbred lines"
 * (Gaynor et al., 2017)
 */
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
int eigen2(arma::vec& eigval, arma::mat& eigvec, arma::mat X){ // Must pass eigval and eigvec by reference or modifications occur on local copy
  char JOBZ = 'V';
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

/*
 * Reads a text file into an arma::mat
 * Must be supplied with the correct number of rows and columns for the matrix
 * A header can be skipped by setting skipRows to 1
 * Row names can be skipped by setting skipCols to 1
 */
arma::mat readMat(std::string fileName, int rows, int cols, char sep=' ',
                  int skipRows=0, int skipCols=0){
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
arma::mat makeX(arma::vec& x){
  int nTrain = x.n_elem;
  double nLevels = x.max();
  arma::mat X(nTrain,nLevels);
  if(nLevels==1){
    X.ones();
  }else{
    X.zeros();
    X.col(0).ones();
    Rcpp::IntegerVector tmp = Rcpp::Range(1,nLevels-1);
    arma::uvec colNeg = Rcpp::as<arma::uvec>(tmp);
    for(int i=0; i<nTrain; ++i){
      if(x(i)==nLevels){
        arma::uvec tmp(1);
        tmp(0) = i;
        X(tmp,colNeg).fill(-1.0);
      }else{
        X(i,int(x(i)))=1;
      }
    }
  }
  return X;
}

// Produces a genotype design matrix
arma::mat makeZ(arma::vec& z, int nPred){
  int nTrain = z.n_elem;
  int nLevels = z.max();
  arma::mat Z(nTrain,nLevels+nPred);
  Z.zeros();
  for(int i=0; i<nTrain; ++i){
    Z(i,int(z(i)-1)) = 1;
  }
  return Z;
}

// Generates weighted matrix
// Allows for heterogenous variance due to unequal replication
arma::mat sweepReps(arma::mat X, const arma::vec& reps){
  arma::mat output(X.n_rows,X.n_cols);
  arma::vec repsT = sqrt(1/reps);
  for(int i=0; i<X.n_cols; ++i){
    output.col(i) = X.col(i)/repsT;
  }
  return output;
}


// Calculates a distance matrix from a marker matrix
// Uses binomial theorem trick
// Inspired by code from:
// http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
// First described here: 
// http://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
arma::mat fastDist(arma::mat& X){
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

arma::mat fastPairDist(arma::mat& X, arma::mat& Y){
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

// Converts two R data.frames into a single mat
arma::mat df2mat(const Rcpp::DataFrame& DF1, 
                 const Rcpp::DataFrame& DF2){
  int nRows1 = DF1.nrows();
  int nRows2 = DF2.nrows();
  arma::mat X((nRows1+nRows2),DF1.size());
  Rcpp::NumericVector tmp1;
  Rcpp::NumericVector tmp2;
  for(int i=0; i<DF1.size(); ++i){
    tmp1 = DF1[i];
    for(int j=0; j<nRows1; ++j){
      X(j,i) = tmp1[j];
    }
    tmp2 = DF2[i];
    for(int k=0; k<nRows2; ++k){
      X(k+nRows1,i) = tmp2[k];
    }
  }
  return X;
}

// Converts one R data.frame into an mat
arma::mat df2mat1(const Rcpp::DataFrame& DF1){
  int nRows1 = DF1.nrows();
  arma::mat X((nRows1),DF1.size());
  Rcpp::NumericVector tmp1;
  for(int i=0; i<DF1.size(); ++i){
    tmp1 = DF1[i];
    for(int j=0; j<nRows1; ++j){
      X(j,i) = tmp1[j];
    }
  }
  return X;
}

// Gaussian kernel function
// D is an Euclidean distance matrix
// theta is the tuning parameter
arma::mat calcGK(arma::mat& D, double theta){
  return exp(-1.0*square(D/theta));
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

// Solve univariate mixed model
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

// Solve multivariate mixed model
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
  bool converging=true;
  while(converging){
    arma::mat VeNew = arma::zeros<arma::mat>(m,m);
    arma::mat VuNew = arma::zeros<arma::mat>(m,m);
    for(int i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*arma::inv_sympd(eigval(i)*Vu+Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    arma::mat BNew = (Yt - Gt)*W;
    arma::mat sigma(m,m);
    for(int i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*arma::inv_sympd(eigval(i)*Vu+Ve+tol*arma::eye(m,m))*(eigval(i)*Vu);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*((Yt.col(i)-BNew*Xt.col(i)-Gt.col(i))*(Yt.col(i)-BNew*Xt.col(i)-Gt.col(i)).t()+sigma);
    }
    double denom = fabs(sum(Ve.diag()));
    if(denom>0.0){
      double numer = fabs(sum(VeNew.diag()-Ve.diag()));
      if((numer/denom)<tol) converging=false;
    }
    Ve = VeNew;
    Vu = VuNew;
    B = BNew;
  }
  arma::mat HI = inv_sympd(kron(ZKZ, Vu)+kron(arma::eye(n,n), Ve)+tol*arma::eye(n*m,n*m));
  arma::mat E = Y.t() - B*X.t();
  arma::mat U = kron(K, Vu)*kron(Z.t(), arma::eye(m,m))*(HI*arma::vectorise(E)); //BLUPs
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

// Objective function for Gaussian kernel method
Rcpp::List gaussObj(double theta, Rcpp::List args){
  Rcpp::List output;
  arma::mat D = args["D"];
  output = solveUVM(args["y"],args["X"],args["Z"],calcGK(D,theta));
  return Rcpp::List::create(Rcpp::Named("objective")=output["LL"],
                            Rcpp::Named("output")=output);
}

// Called by GK_AS R function
// [[Rcpp::export(.callGK_AS)]]
Rcpp::List callGK_AS(arma::mat y, arma::vec x, arma::vec reps,
                     std::string genoTrain, int nMarker, double maxTheta, 
                     int maxIter, bool writeForPred=true){
  Rcpp::List output;
  arma::mat X = makeX(x);
  arma::mat Z;
  Z.eye(y.n_elem,y.n_elem);
  y = sweepReps(y,reps);
  X = sweepReps(X,reps);
  Z = sweepReps(Z,reps);
  arma::mat M = readMat(genoTrain,y.n_elem,nMarker,',',0,1);
  arma::mat D = fastDist(M);
  output = optimize(*gaussObj,Rcpp::List::create(Rcpp::Named("y")=y,
                                                 Rcpp::Named("X")=X,
                                                 Rcpp::Named("Z")=Z,
                                                 Rcpp::Named("D")=D),
                                                 1e-10,
                                                 maxTheta*D.max(),
                                                 maxIter,
                                                 true);
  if(writeForPred){
    M.save("M.bin");
    arma::mat theta(1,1);
    theta(0,0) = double(output["parameter"]);
    theta.save("theta.bin");
    Rcpp::List tmp = output["output"];
    double Vu = tmp["Vu"];
    double Ve = tmp["Ve"];
    Rcpp::NumericMatrix tmpBeta = tmp["beta"];
    arma::mat beta = Rcpp::as<arma::mat>(tmpBeta);
    arma::mat intercept(1,1);
    intercept = beta(0,0);
    intercept.save("intercept.bin");
    arma::mat Yt = y - X*beta;
    arma::mat W = inv(calcGK(D, double(theta(0,0)))+Ve/Vu*arma::eye(D.n_rows,D.n_cols));
    W = W*Yt;
    W.save("W.bin");
  }
  return output;
}

// Predicts BLUP from callGK_AS output
// Called from PredGK_AS R function
// [[Rcpp::export(.callPredGK_AS)]]
arma::mat callPredGK_AS(const Rcpp::DataFrame& genoPred){
  arma::mat X;
  X = df2mat1(genoPred);
  arma::mat M;
  M.load("M.bin");
  arma::mat W;
  W.load("W.bin");
  arma::mat tmpInter;
  tmpInter.load("intercept.bin");
  double intercept = tmpInter(0,0);
  arma::mat tmpTheta;
  tmpTheta.load("theta.bin");
  double theta = tmpTheta(0,0);
  arma::mat Cut = fastPairDist(X,M);
  Cut = calcGK(Cut,theta);
  return Cut*W+intercept;
}

// Called from RRBLUP_AS R function
// [[Rcpp::export(.callRRBLUP_AS)]]
Rcpp::List callRRBLUP_AS(arma::mat y, arma::vec x, arma::vec reps, 
                         std::string genoTrain, int nMarker){
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,y.n_elem,nMarker,',',0,1);
  y = sweepReps(y,reps);
  X = sweepReps(X,reps);
  Z = sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  return solveUVM(y, X, Z, K);
}

// Called from RRBLUP_AS_MV R function
// [[Rcpp::export(.callRRBLUP_AS_MV)]]
Rcpp::List callRRBLUP_AS_MV(arma::mat Y, arma::vec x, arma::vec reps, 
                            std::string genoTrain, int nMarker){
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,Y.n_rows,nMarker,',',0,1);
  Y = sweepReps(Y,reps);
  X = sweepReps(X,reps);
  Z = sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  return solveMVM(Y, X, Z, K);
}

// Called from QUICK_AS R function
// [[Rcpp::export(.callQUICK_AS)]]
Rcpp::List callQUICK_AS(arma::mat y, arma::vec x, arma::vec reps, 
                        std::string genoTrain, int nMarker, double varA,
                        double varE){
  int n=y.n_rows;
  double delta=varE/(varA/double(nMarker));
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,y.n_elem,nMarker,',',0,1);
  y = sweepReps(y,reps);
  X = sweepReps(X,reps);
  Z = sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  arma::mat ZK = Z*K;
  arma::mat Hinv = inv(ZK*Z.t()+delta*arma::eye(n,n));
  arma::mat XHinv = X.t()*Hinv;
  arma::mat beta = solve(XHinv*X,XHinv*y);
  arma::mat u = ZK.t()*(Hinv*(y-X*beta));
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u);
}
