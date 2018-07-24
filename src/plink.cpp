// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>

// Write PLINK PED file
// [[Rcpp::export]]
void writePlinkPed(Rcpp::DataFrame & fam,
                   const arma::field<arma::Cube<unsigned char> > & geno,
                   const arma::uvec & lociPerChr,
                   arma::uvec & lociLoc,
                   const std::string file) {
  
  // R to C++ index correction
  lociLoc -= 1;

  int nInd   = geno(0).n_slices;
  int nChr   = geno.n_elem;
  
  Rcpp::IntegerVector family = fam["family"];
  Rcpp::IntegerVector id     = fam["id"];
  Rcpp::IntegerVector father = fam["father"];
  Rcpp::IntegerVector mother = fam["mother"];
  Rcpp::IntegerVector gender = fam["gender"];
  Rcpp::NumericVector pheno  = fam["pheno"];

  std::ofstream plinkPed;
  plinkPed.open(file, std::ofstream::trunc);
  for (int ind = 0; ind < nInd; ind++) {
    // Family/Individual data
    plinkPed << family[ind] << " " << id[ind]     << " " << father[ind] << " "
             << mother[ind] << " " << gender[ind] << " " << pheno[ind];
    // Genotype data
    int loc = 0;
    for (int chr = 0; chr < nChr; ++chr) {
      if (lociPerChr(chr) > 0) {
        for (int locInChr = 1; locInChr <= lociPerChr(chr); ++locInChr) {
          /*
            std::cout << "Ind "       << ind      << " Chr "      << chr          << " LocOverall " << loc;
            std::cout << " LocOnChr " << locInChr << " PosOnChr " << lociLoc(loc) << " Alleles ";
            std::cout << geno(chr).slice(ind).col(0).row(lociLoc(loc)) << " "
                      << geno(chr).slice(ind).col(1).row(lociLoc(loc)) << " " << std::endl;
          */
          // arma::conv_to<int>::from() needed due to file formatting
          // @todo can we make this more efficient by not calling the conversion function so many times?
          plinkPed << " " << arma::conv_to<int>::from(geno(chr).slice(ind).col(0).row(lociLoc(loc)) + 1)
                   << " " << arma::conv_to<int>::from(geno(chr).slice(ind).col(1).row(lociLoc(loc)) + 1);
          loc += 1;
        }
      }
      
    }
    plinkPed << std::endl;
  }
  plinkPed.close();
}

/*** R
writePlinkPed(fam = fam,
              geno = pop@geno,
              lociPerChr = tmp$lociPerChr,
              lociLoc = tmp$lociLoc,
              file = paste0(baseName, ".ped"))
*/