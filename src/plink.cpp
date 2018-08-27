// [[Rcpp::depends(RcppArmadillo)]]
#include "alphasimr.h"

// Write PLINK PED file
// [[Rcpp::export]]
void writePlinkPed(const Rcpp::DataFrame                         & fam,
                   const arma::field<arma::Cube<unsigned char> > & geno,
                   const arma::uvec                              & lociPerChr,
                         arma::uvec                                lociLoc,
                   const Rcpp::String                            & file) {
  
  // R to C++ index correction
  lociLoc -= 1;

  int nInd = geno(0).n_slices;
  int nChr = geno.n_elem;
  
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
        // geno stores 0 and 1 as "unsigned char" [0 - 255]
        // to get ASCII characters 1 and 2 we need values 49 and 50 http://www.asciitable.com
        // and cast these "unsigned char" values as "char"
        arma::Col<unsigned char> gamete1 = geno(chr).slice(ind).col(0) + 49;
        arma::Col<unsigned char> gamete2 = geno(chr).slice(ind).col(1) + 49;
        for (int locInChr = 1; locInChr <= lociPerChr(chr); ++locInChr) {
          /*
            std::cout << "Ind "       << ind      << " Chr "      << chr          << " LocOverall " << loc
                      << " LocOnChr " << locInChr << " PosOnChr " << lociLoc(loc) << " Alleles "
                      << " " << char(gamete1(lociLoc(loc)))
                      << " " << char(gamete2(lociLoc(loc))) << "\n";
          */
          plinkPed << " " << char(gamete1(lociLoc(loc)))
                   << " " << char(gamete2(lociLoc(loc)));
          loc += 1;
        }
      }
    }
    plinkPed << "\n";
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
