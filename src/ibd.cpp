#include "alphasimr.h"

// Converts founder IBD into hap format for SimParam
// [[Rcpp::export]]
arma::field< //individual
  arma::field< //chromosome
    arma::field< //ploidy
      arma::Mat<int> > > > getFounderIbd(const arma::field<arma::ivec>& founder, 
                                         arma::uword nChr){
        
        // Generate object for output
        arma::field< //individual
          arma::field< //chromosome
            arma::field< //ploidy
              arma::Mat<int> > > > output;
        output.set_size(founder.n_elem);
        
        // Allocate output object
        for(arma::uword i=0; i<founder.n_elem; ++i){
          output(i).set_size(nChr);
          for(arma::uword j=0; j<nChr; ++j){
            output(i)(j).set_size(founder(i).n_elem);
          }
          for(arma::uword k=0; k<founder(i).n_elem; ++k){
            arma::Mat<int> tmp(1,2,arma::fill::ones);
            tmp(0,0) = founder(i)(k);
            for(arma::uword j=0; j<nChr; ++j){
              output(i)(j)(k) = tmp;
            }
          }
        }
        return output;
      }

// Calculates IBD for individual using recombination data and parental IBD
// [[Rcpp::export]]
arma::field< //chromosome
  arma::field< //ploidy
    arma::Mat<int> > > getNonFounderIbd(
        const arma::field<arma::field<arma::Mat<int> > >& recHist,
        const arma::field<arma::field<arma::Mat<int> > >& mother,
        const arma::field<arma::field<arma::Mat<int> > >& father){
      arma::uword nChr = recHist.n_elem;
      arma::uword motherPloidy = mother(0).n_elem;
      arma::uword ploidy = recHist(0).n_elem;
      arma::field<arma::field<arma::Mat<int> > > output;
      output.set_size(nChr);
      
      for(arma::uword i=0; i<nChr; ++i){
        output(i).set_size(ploidy);
        for(arma::uword j=0; j<ploidy; ++j){
          arma::uword nCO = recHist(i)(j).n_rows - 1;
          if(j<motherPloidy/2){
            // Pull from mother
            if(nCO==0){
              // Direct copy of chromosome IBD
              output(i)(j) = mother(i)(recHist(i)(j)(0,0)-1);
            }else{
              arma::Mat<int> X;
              // Transfer first section
              X = mother(i)(recHist(i)(j)(0,0)-1);
              X = X.rows(find(X.col(1)<recHist(i)(j)(1,1)));
              
              // Transfer middle section(s)?
              if(nCO>1){
                for(arma::uword k=1; k<nCO; ++k){
                  arma::Mat<int> Y = mother(i)(recHist(i)(j)(k,0)-1);
                  // Trim end
                  Y = Y.rows(find(Y.col(1)<recHist(i)(j)(k+1,1)));
                  
                  // Trim front
                  arma::uword l;
                  for(l=0; l<Y.n_rows; l++){
                    if(Y(l,1)>recHist(i)(j)(k+1,1)){
                      --l;
                      break;
                    }
                  }
                  Y = Y.rows(l,Y.n_rows-1);
                  
                  // Check for haplotype match
                  if(Y(0,1)==X(X.n_rows-1,1)){
                    // Drop first row
                    Y.shed_row(0);
                  }else{
                    // Change first site to crossover location
                    Y(0,1) = recHist(i)(j)(k,1);
                  }
                  
                  X = join_cols(X, Y);
                }
              }
              
              // Transfer last section
              arma::Mat<int> Y = mother(i)(recHist(i)(j)(nCO,0)-1);
              
              // Trim front
              arma::uword l;
              arma::uword m = Y.n_rows-1;
              for(l=0; l<Y.n_rows; l++){
                if(Y(l,1)>recHist(i)(j)(nCO,1)){
                  m = l-1;
                  break;
                }
              }
              Y = Y.rows(m,Y.n_rows-1);
              
              // Check for haplotype match
              if(Y(0,1)==X(X.n_rows-1,1)){
                // Drop first row
                Y.shed_row(0);
              }else{
                // Change first site to crossover location
                Y(0,1) = recHist(i)(j)(nCO,1);
              }
              
              output(i)(j) = join_cols(X, Y);
            }
          }else{
            // Pull from father
            if(nCO==0){
              // Direct copy of chromosome IBD
              output(i)(j) = father(i)(recHist(i)(j)(0,0)-1 );
            }else{
              arma::Mat<int> X;
              // Transfer first section
              X = father(i)(recHist(i)(j)(0,0)-1);
              X = X.rows(find(X.col(1)<recHist(i)(j)(1,1)));
              
              // Transfer middle section(s)?
              if(nCO>1){
                for(arma::uword k=1; k<nCO; ++k){
                  arma::Mat<int> Y = father(i)(recHist(i)(j)(k,0)-1);
                  // Trim end
                  Y = Y.rows(find(Y.col(1)<recHist(i)(j)(k+1,1)));
                  
                  // Trim front
                  arma::uword l;
                  for(l=0; l<Y.n_rows; l++){
                    if(Y(l,1)>recHist(i)(j)(k+1,1)){
                      --l;
                      break;
                    }
                  }
                  Y = Y.rows(l,Y.n_rows-1);
                  
                  // Check for haplotype match
                  if(Y(0,1)==X(X.n_rows-1,1)){
                    // Drop first row
                    Y.shed_row(0);
                  }else{
                    // Change first site to crossover location
                    Y(0,1) = recHist(i)(j)(k,1);
                  }
                  
                  X = join_cols(X, Y);
                }
              }
              
              // Transfer last section
              arma::Mat<int> Y = father(i)(recHist(i)(j)(nCO,0)-1);
              
              // Trim front
              arma::uword l;
              arma::uword m = Y.n_rows-1;
              for(l=0; l<Y.n_rows; l++){
                if(Y(l,1)>recHist(i)(j)(nCO,1)){
                  m = l-1;
                  break;
                }
              }
              Y = Y.rows(m,Y.n_rows-1);
              
              // Check for haplotype match
              if(Y(0,1)==X(X.n_rows-1,1)){
                // Drop first row
                Y.shed_row(0);
              }else{
                // Change first site to crossover location
                Y(0,1) = recHist(i)(j)(nCO,1);
              }
              
              output(i)(j) = join_cols(X, Y);
            }
          }
        }
      }
      
      return output;
    }
