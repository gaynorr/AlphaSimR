#ifndef OPTIMIZE_H
#define OPTIMIZE_H

Rcpp::List optimize(Rcpp::List (*objective)(double, Rcpp::List), 
                    Rcpp::List args, double l, double u, int maxIter=1000, 
                    bool maximize=false, bool evalU = false, bool evalL = false, 
                    double eps=1.0e-9);

#endif