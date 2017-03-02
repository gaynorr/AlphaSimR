// A stop-gap solution until MaCS is integrated

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @internal
// [[Rcpp::export]]
arma::Mat<unsigned char> readAF(int nHap, int segSites, arma::ivec keep){
  keep -= 1; //R to C++
  arma::Mat<unsigned char> output(nHap, keep.n_elem);
  char sep = ' ';
  std::string iFileName = "MacsHaplotypes.txt";
  std::ifstream iFile(iFileName.c_str());
  std::string line;
  //Read rows
  for(int i=0; i<nHap; ++i){
    std::getline(iFile,line);
    std::stringstream lineStream(line);
    std::string cell;
    //Skip blank column
    std::getline(lineStream,cell,sep);
    //Read columns
    int k=0;
    for(int j=0; j<segSites; ++j){
      std::getline(lineStream,cell,sep);
      if(j==keep[k]){
        output(i,k) = std::atof(cell.c_str());
        ++k;
      }
    }
  }
  iFile.close();
  return output;
}
