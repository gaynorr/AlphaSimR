#include "alphasimr.h"

// Write PLINK PED file
// [[Rcpp::export]]
void writePlinkPed(const Rcpp::DataFrame          & fam,
                   const arma::Mat<unsigned char> & haplo,
                   const arma::uword              & nInd,
                   const arma::uword              & ploidy,
                   const arma::uword              & nLoc,
                   const Rcpp::String             & file) {
  
  // FAM file stuff  
  Rcpp::IntegerVector family = fam["family"];
  Rcpp::IntegerVector id     = fam["id"];
  Rcpp::IntegerVector father = fam["father"];
  Rcpp::IntegerVector mother = fam["mother"];
  Rcpp::IntegerVector sex = fam["sex"];
  Rcpp::NumericVector pheno  = fam["pheno"];
  
  std::ofstream plinkPed;
  plinkPed.open(file, std::ofstream::trunc);
  arma::Col<char> alleles(nLoc * ploidy);
  for (arma::uword ind = 0; ind < nInd; ind++) {
    // Family & Individual data
    plinkPed << family[ind] << " " << id[ind]     << " " << father[ind] << " "
             << mother[ind] << " " << sex[ind] << " " << pheno[ind];
    // Haplotype data
    // ... first gather locus alleles from all chromosomes
    arma::uword indLastChrom = (ind + 1) * ploidy;
    for (arma::uword p = ploidy; p > 0; p--) {
      arma::uword indChrom = indLastChrom - p;
      for (arma::uword loc = 0; loc < nLoc; loc++) {
        arma::uword locLastChrom = (loc + 1) * ploidy;
        /*
         std::cout << " Ind "       << ind + 1
                   << " Chrom "     << ploidy - p + 1
                   << " Locus "     << loc + 1
                   << " Haplo row " << indChrom
                   << " Allele "    << haplo(indChrom, loc) + 49
                   << " Allele "    << char(haplo(indChrom, loc) + 49)
                   << "\n";
         */
        // haplo stores 0 and 1 as "unsigned char" [0 - 255]
        // to get ASCII characters 1 and 2 we need values 49 and 50 http://www.asciitable.com
        // and cast these "unsigned char" values as "char"
        alleles(locLastChrom - p) = char(haplo(indChrom, loc) + 49);
      }
    }
    // ... now print locus alleles from all chromosomes
    for (arma::uword loc = 0; loc < (nLoc * ploidy); loc++) {
      plinkPed << " " << alleles(loc);
    }
    plinkPed << "\n";
  }
  plinkPed.close();
}
