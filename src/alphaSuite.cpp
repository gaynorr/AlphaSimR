// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"
#include <iostream>
#include <fstream>

// [[Rcpp::export]]
void writeASGenotypes(const arma::Cube<unsigned char> & g,
                      const arma::field<arma::uvec> &locations,
                      const arma::uvec &allLocations,
                      const arma::vec & snpchips, 
                      const std::vector<std::string> & names, 
                      const char missing,
                      const std::string fname){
  
  std::ofstream ASout;
  ASout.open(fname, std::ios::trunc);

  char name[20];
  for (arma::uword i = 0; i < snpchips.n_rows; i++){
    arma::Col<unsigned char> all0 = g.slice(i).col(0);
    arma::Col<unsigned char> all1 = g.slice(i).col(1);
    arma::Col<unsigned char> selected0;
    arma::Col<unsigned char> selected1;
    
    if (snpchips(i) == 0) {
      selected0 = all0.elem(allLocations - 1);
      selected1 = all1.elem(allLocations - 1);      
    }
    else{
      selected0 = all0.elem(locations(snpchips(i) - 1) - 1);
      selected1 = all1.elem(locations(snpchips(i) - 1) - 1);
    }
    
    arma::Col<unsigned char> selectedg = selected0 + selected1;
    
    std::sprintf(name,"%s",names[i].c_str());
    ASout << name;
    
    if (snpchips(i) == 0) {
      for (arma::uword k = 0; k < allLocations.n_rows; k++) {
        ASout << " " << char(selectedg(k) + 48);
      }      
    }
    else
    {
      arma::uword cur = 0;
      for (arma::uword k = 0; k < allLocations.n_rows; k++) {
        if ( (cur < locations(snpchips(i) - 1).n_rows) && ((allLocations(k) - 1) == (locations(snpchips(i) - 1)(cur) - 1)) ){
          ASout << " " << char(selectedg(cur) + 48);
          cur ++;
        }
        else {
          ASout << " " << missing;
        }
      }
    }
    
    ASout << "\n";
  }
  
  ASout.close();
}

// [[Rcpp::export]]
void writeASHaplotypes(const arma::Cube<unsigned char> & g,
                      const arma::field<arma::uvec> &locations,
                      const arma::uvec &allLocations,
                      const arma::vec & snpchips, 
                      const std::vector<std::string> & names,
                      const char missing,
                      const std::string fname){
  
  std::ofstream ASout;
  ASout.open(fname, std::ios::trunc);
  
  char name[20];
  for (arma::uword i = 0; i < snpchips.n_rows; i++){
    for (arma::uword j = 0; j < 2; j ++){
      arma::Col<unsigned char> all = g.slice(i).col(j);

      std::sprintf(name,"%s",names[i].c_str());
      ASout << name;

      if (snpchips(i) == 0) {
        arma::Col<unsigned char> selected = all.elem(allLocations - 1);
        for (arma::uword k = 0; k < allLocations.n_rows; k++) {
          ASout << " " << char(selected(k) + 48);
        }      
      }
      else {
        arma::uword cur = 0;
        for (arma::uword k = 0; k < allLocations.n_rows; k++) {
          arma::Col<unsigned char> selected = all.elem(locations(snpchips(i) - 1) - 1);
          if ( (cur < locations(snpchips(i) - 1).n_rows) && ((allLocations(k) - 1) == (locations(snpchips(i) - 1)(cur) - 1)) ){
            ASout << " " << char(selected(cur) + 48);
            cur ++;
          }
          else {
            ASout << " " << missing;
          }
        }
      }
      
      ASout << "\n";
    }
  }
  
  ASout.close();
}