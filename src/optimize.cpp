#include <Rcpp.h>
using namespace Rcpp;
/*
 * Modified version of Brent's method for numerical optimization
 * objective is a pointer to a function returning a List
 *   objective takes a double and a List as arguments
 *   optimization is for the double
 *   the returned list includes "objective" and list of additional output
 * args is a list which is passed to the objective
 * l is the lower bound of optimization
 * u is the upper bound of optimization
 * maxIter is the maximum number of iterations
 *   the best obtained value is returned when maxIter is reached
 * if firstLast is true, l and u are evaluated if convergence isn't reached
 * maximize searches for a maximum value
 * eps controls the tolerance for convergence
 */
List optimize(List (*objective)(double, List), List args, double l, 
              double u, int maxIter=1000, bool maximize=false, 
              bool evalU=false, bool evalL=false, double eps=1.0e-9){
  double MACHEPS_SQRT = sqrt(2.2204460492503131e-016);
  double c = (3.0-std::sqrt(double(5.0)))/2.0;
  double x = l+c*(u-l);
  double v = x;
  double w = x;
  double e = 0.0;
  double lInt = l;
  double uInt = u;
  List fOut = objective(x, args);
  List output = fOut["output"];
  double fx = fOut["objective"];
  if(maximize) fx = -fx;
  double fv = fx;
  double fw = fx;
  
  int numiter = 0;
  
  while(true){
    double m = 0.5*(l+u);
    double tol = MACHEPS_SQRT*std::abs(x)+eps;
    double tol2 = 2.0*tol;
    
    // Check the stopping criterion
    if(std::abs(x-m)<=(tol2-(0.5*(u-l)))){
      break;
    }
    
    // Check maximum iterations
    if (++numiter>maxIter){
      break;
    }
    
    double p=0,q=0,r=0,d=0,z=0;
    
    if(std::abs(e)>tol){
      // Fit parabola
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if(q>0.0){
        p = -p;
      }else{
        q = -q;
      }
      r = e;
      e = d;
    }
    
    if((std::abs(p)<std::abs(0.5*q*r))&&
       (p<(q*(l-x)))&&
       (p<(q*(u-x)))){
      // Parabolic interpolation step
      d = p/q;
      z = x+d;
      // objective must not be evaluated too close to l or u
      if( ((z-l)<tol2)||((u-z)<tol2) ){
        d = (x<m)?tol:-tol;
      }
    }else{
      // Golden section step
      e = (x<m)?(u-x):(l-x);
      d = c*e;
    }
    
    // objective must not be evaluated too close to x
    if(std::abs(d)>=tol){
      z = x+d;
    }else if(d>0.0){
      z = x+tol;
    }else{
      z = x-tol;
    }
    fOut = objective(z,args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    
    // Update
    if(funcu<=fx){
      if(z<x){
        u = x;
      }else{
        l = x;
      }
      
      v = w; 
      fv = fw;
      w = x; 
      fw = fx;
      x = z;
      output = fOut["output"];
      fx = funcu;
    }else{
      if(z<x){
        l = z;
      }else{
        u = z;
      }
      
      if( (funcu<=fw)||(w==x) ){
        v = w;
        fv = fw;
        w = z;
        fw = funcu;
      }else if( (funcu<=fv)||(v==x)||(v==w) ){
        v = z;
        fv = funcu;
      }
    }
  }
  if(evalU){
    fOut = objective(uInt,args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    if(funcu<fx){
      fx = funcu;
      output = fOut["output"];
      x = uInt;
    }
  }
  if(evalL){
    fOut = objective(lInt,args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    if(funcu<fx){
      fx = funcu;
      output = fOut["output"];
      x = lInt;
    }
  }
  if(maximize) fx = -fx;
  return List::create(Named("parameter") = x,
                      Named("objective") = fx,
                      Named("output") = output);
}
