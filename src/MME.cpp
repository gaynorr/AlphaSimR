/*
 * Initial code from:
 * "A two-part strategy for using genomic selection to develop inbred lines"
 * (Gaynor et al., 2017)
 * 
 * Some of functions are not used
 * 
 * Based on code in R/EMMREML and R/rrBLUP
 * 
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

// A class used by solveMKM
// This class is used to avoid copying data to R when calling optim
class DataMKM{
public:
  arma::mat y;
  arma::mat X;
  arma::field<arma::mat> Zlist;
  arma::field<arma::mat> Klist;
  int n;
  int q;
  int nre; //Number of random effects
  double df;
  double offset;
  arma::mat ZK;
  arma::mat S;
  arma::mat ZKZ;
  arma::mat SZKZ;
  arma::vec weights;
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec eta;
  double delta;
  double ll;
  arma::mat u;
  arma::mat beta;
  double Vu;
  double Ve;
  
  // Constructor
  DataMKM(arma::mat& y_, arma::mat& X_, 
          arma::field<arma::mat>& Zlist_,
          arma::field<arma::mat>& Klist_) :
    y(y_), X(X_), Zlist(Zlist_), Klist(Klist_) {
    n = y_.n_rows;
    q = X_.n_cols;
    nre = Zlist_.n_elem;
    df = double(n)-double(q);
    offset = log(double(n));
    S.set_size(n,n);
    S = arma::eye(n,n) - X*arma::inv_sympd(X.t()*X)*X.t();
    ZKZ.set_size(n,n);
    SZKZ.set_size(n,n);
    eta.set_size(n);
  }
  
  // Sets ZKZ matrix
  void setZKZ(Rcpp::NumericVector x){
    weights = Rcpp::as<arma::vec>(x);
    weights = weights/sum(weights);
    ZKZ = weights(0)*Zlist(0)*Klist(0)*Zlist(0).t();
    for(int i=1; i<nre; ++i){
      ZKZ += weights(i)*Zlist(i)*Klist(i)*Zlist(i).t();
    }
    SZKZ = S*(ZKZ+offset*arma::eye(n,n))*S;
  }
  
  // Decomposes SZKZ matrix
  Rcpp::NumericVector decompSZKZ(){
    eigval.set_size(n);
    eigvec.set_size(n,n);
    eigen2(eigval, eigvec, SZKZ);
    // Drop eigenvalues
    eigval = eigval(arma::span(q,eigvec.n_cols-1)) - offset;
    eigvec = eigvec(arma::span(0,eigvec.n_rows-1),
                    arma::span(q,eigvec.n_cols-1));
    eta = eigvec.t()*y;
    // Calculate variance
    Rcpp::List optRes = optimize(*objREML,
                                 Rcpp::List::create(
                                   Rcpp::Named("df")=df,
                                   Rcpp::Named("eta")=eta,
                                   Rcpp::Named("lambda")=eigval), 
                                   1.0e-10, 1.0e10);
    delta = optRes["parameter"];
    ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
    Rcpp::NumericVector output = optRes["objective"];
    return output;
  }
  
  void setZK(){
    int nColZ=0;
    int nRowK=0;
    int nColK=0;
    for(int i=0; i<nre; ++i){
      nColZ += Zlist(i).n_cols;
      nRowK += Klist(i).n_rows;
      nColK += Klist(i).n_cols;
    }
    arma::mat Z(n,nColZ);
    arma::mat K(nRowK,nColK,arma::fill::zeros);
    nColZ = 0;
    nRowK = 0;
    nColK = 0;
    for(int i=0; i<nre; ++i){
      Z(arma::span(0,n-1),arma::span(nColK,Zlist(i).n_cols-1+nColK)) = Zlist(i);
      nColZ += Zlist(i).n_cols;
      K(arma::span(nRowK,Klist(i).n_rows-1+nRowK),arma::span(nColK,Klist(i).n_cols-1+nColK)) = weights(i)*Klist(i);
      nColK += Klist(i).n_cols;
      nRowK += Klist(i).n_rows;
    }
    ZK = Z*K;
  }
  
  void solveMME(){
    arma::mat Hinv = arma::inv_sympd(ZKZ+delta*arma::eye(n,n));
    arma::mat XHinv = X.t()*Hinv;
    beta = solve(XHinv*X,XHinv*y);
    u = ZK.t()*(Hinv*(y-X*beta));
    Vu = sum(eta%eta/(eigval+delta))/df;
    Ve = delta*Vu;
  }
};

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
    //Assuming there won't be any Inf values
    D.elem(find_nonfinite(D)).fill(0.0); 
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
    //Assuming there won't be any Inf values
    D.elem(find_nonfinite(D)).fill(0.0); 
  }
  return D; 
}

// Calculates VanRaden's G matrix
// uses 2,1,0 coding of markers
// modifies X and p
arma::mat calcG(arma::mat& X, arma::rowvec& p){
  p = mean(X,0)/2.0;
  X.each_row() -= 2*p;
  arma::mat G = X*X.t();
  G = G/(2.0*sum(p%(1-p)));
  return G;
}

// Creates G matrix from markers using readMat for markers
arma::mat makeG(std::string fileName, int rows, int cols, char sep=' ',
                int skipRows=0, int skipCols=0){
  arma::mat X = readMat(fileName,rows,cols,sep,skipRows,skipCols);
  arma::rowvec p(X.n_cols);
  arma::mat G = calcG(X,p);
  return G;
}

// Gaussian kernel function
// D is an Euclidean distance matrix
// theta is the tuning parameter
arma::mat calcGK(arma::mat& D, double theta){
  return exp(-1.0*square(D/theta));
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
  bool converging=true;
  while(converging){
    arma::mat VeNew = arma::zeros<arma::mat>(m,m);
    arma::mat VuNew = arma::zeros<arma::mat>(m,m);
    for(int i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*arma::inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    arma::mat BNew = (Yt - Gt)*W;
    arma::mat sigma(m,m);
    for(int i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*arma::inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(eigval(i)*Vu);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*((Yt.col(i)-BNew*Xt.col(i)-Gt.col(i))*
        (Yt.col(i)-BNew*Xt.col(i)-Gt.col(i)).t()+sigma);
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

// Objective function for calculating variance component weights
// Used by the R function optim
// [[Rcpp::export]]
Rcpp::NumericVector objWeights(Rcpp::NumericVector x, 
                               SEXP ptrData){
  Rcpp::XPtr<DataMKM> ptr = Rcpp::as<Rcpp::XPtr<DataMKM> >(ptrData);
  ptr->setZKZ(x);
  Rcpp::NumericVector output = ptr->decompSZKZ();
  return output;
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
  Rcpp::XPtr<DataMKM> ptrData(new DataMKM(y, X, Zlist,Klist), 
                              true);
  //Calculate weights
  Rcpp::Function optim("optim");
  Rcpp::Environment E("package:AlphaSimR");
  Rcpp::Function objWeightsR = E["objWeightsR"];
  Rcpp::NumericVector init(ptrData->nre,1.0/double(ptrData->nre));
  Rcpp::NumericVector zeros(ptrData->nre); // Lower bounds
  Rcpp::NumericVector ones(ptrData->nre,1.0); // Upper bounds
  optim(Rcpp::_["par"]=init, Rcpp::_["fn"]=objWeightsR, 
        Rcpp::_["ptrData"]=ptrData, Rcpp::_["method"]="L-BFGS-B", 
        Rcpp::_["lower"]=zeros, Rcpp::_["upper"]=ones);
  
  //Find solution
  ptrData->setZK();
  ptrData->solveMME();
  return Rcpp::List::create(Rcpp::Named("Vu")=ptrData->Vu,
                            Rcpp::Named("Ve")=ptrData->Ve,
                            Rcpp::Named("beta")=ptrData->beta,
                            Rcpp::Named("u")=ptrData->u,
                            Rcpp::Named("LL")=ptrData->ll,
                            Rcpp::Named("weights")=ptrData->weights);
}

// Objective function for Gaussian kernel method
Rcpp::List gaussObj(double theta, Rcpp::List args){
  Rcpp::List output;
  arma::mat D = args["D"];
  output = solveUVM(args["y"],args["X"],args["Z"],calcGK(D,theta));
  return Rcpp::List::create(Rcpp::Named("objective")=output["LL"],
                            Rcpp::Named("output")=output);
}

// Objective function for multivariate Gaussian kernel method
Rcpp::List gaussObjMV(double theta, Rcpp::List args){
  Rcpp::List output;
  arma::mat D = args["D"];
  output = solveMVM(args["Y"],args["X"],args["Z"],calcGK(D,theta));
  return Rcpp::List::create(Rcpp::Named("objective")=output["LL"],
                            Rcpp::Named("output")=output);
}

// Called by GK R function
// [[Rcpp::export]]
Rcpp::List callGK(arma::mat y, arma::uvec x, arma::vec reps,
                     std::string genoTrain, int nMarker, 
                     double maxTheta, int maxIter, 
                     bool writeForPred=true){
  Rcpp::List output;
  arma::mat X = makeX(x);
  arma::mat Z;
  Z.eye(y.n_elem,y.n_elem);
  reps = sqrt(1.0/reps);
  sweepReps(y,reps);
  sweepReps(X,reps);
  sweepReps(Z,reps);
  arma::mat M = readMat(genoTrain,y.n_elem,nMarker,' ');
  arma::mat D = fastDist(M);
  output = optimize(*gaussObj,
                    Rcpp::List::create(Rcpp::Named("y")=y,
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
    arma::mat W = inv(calcGK(D, double(theta(0,0)))+
      Ve/Vu*arma::eye(D.n_rows,D.n_cols));
    W = W*Yt;
    W.save("W.bin");
  }
  return output;
}

// Predicts BLUP from callGK output
// Called from PredGK R function
// [[Rcpp::export]]
arma::mat callPredGK(const Rcpp::DataFrame& genoPred){
  arma::mat X;
  //X = df2mat1(genoPred);
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

// Called by RRBLUP function
// [[Rcpp::export]]
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x, arma::vec reps, 
                         std::string genoTrain, int nMarker){
  arma::mat X = makeX(x);
  arma::mat Z = readMat(genoTrain,y.n_elem,nMarker,' ');
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
  arma::mat Z = readMat(genoTrain,Y.n_rows,nMarker,' ');
  reps = sqrt(1.0/reps);
  sweepReps(Y,reps);
  sweepReps(X,reps);
  sweepReps(Z,reps);
  arma::mat K = arma::eye(nMarker,nMarker);
  return solveMVM(Y, X, Z, K);
}

Rcpp::List callSC_GBLUP(arma::mat y, arma::uvec x, arma::vec reps,
                        std::string genoFemale, int nFemale, 
                        arma::uvec fPar, std::string genoMale, 
                        int nMale, arma::uvec mPar, int nMarker){
  int n = y.n_rows;
  arma::mat X = makeX(x);
  arma::field<arma::mat> Zlist(3);
  Zlist(0) = makeZ(fPar, nFemale);
  Zlist(1) = makeZ(mPar, nMale);
  Zlist(2) = arma::eye(n,n);
  arma::field<arma::mat> Klist(3);
  Klist(0) = makeG(genoFemale,nFemale,nMarker);
  Klist(1) = makeG(genoMale,nMale,nMarker);
  // K3 = Z1*K1*Z1' % Z2*K2*Z2'
  Klist(2) = (Zlist(0)*Klist(0)*Zlist(0).t())%
  (Zlist(1)*Klist(1)*Zlist(1).t());
  reps = sqrt(1.0/reps);
  sweepReps(y, reps);
  sweepReps(X, reps);
  sweepReps(Zlist(0), reps);
  sweepReps(Zlist(1), reps);
  sweepReps(Zlist(2), reps);
  return solveMKM(y,X,Zlist,Klist);
}

