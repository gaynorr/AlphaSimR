#include "alphasimr.h"
#include "pedigree.h"
#include <algorithm>
#include <cmath>
#include <utility>

// Numerator relationships from a pedigree: inbreeding coefficients and
// blocks of the relationship matrix by the indirect method of Colleau
// (2002, doi:10.1186/1297-9686-34-4-409), the inverse of the relationship
// matrix by the rules of Henderson (1976,
// doi:10.2307/2529339) and Quaas (1976, doi:10.2307/2529279), and the
// relationship matrix itself by Henderson's tabular method.
//
// For a ploidy level above 2, inheritance is taken to be autopolyploid
// without double reduction, following Kerr et al. (2012,
// doi:10.1007/s00122-012-1785-y) and Hamilton and Kerr (2018,
// doi:10.1007/s00122-017-3041-y), with relationships scaled so that an
// individual's diagonal is 1+(v-1)F for ploidy v. Off-diagonals are then
// the average of the parents' at every ploidy, so only the inbreeding
// recursion depends on it.
//
// All three are written for the decomposition A = T*diag(D)*T', in which
// each individual's breeding value is the average of its parents' plus a
// Mendelian sampling term of variance D. A doubled haploid fits this with
// its one parent recorded as both mother and father: the two halves of the
// parent's share add up to the whole of it, which is right, because a
// doubled haploid's genome is two copies of one gamete of its parent. Only
// its Mendelian sampling variance differs from that of a selfed individual,
// so that is the one place where doubled haploids are treated separately.

namespace {

// Mendelian sampling variances at or below this are treated as zero. Exact
// zeros arise from fully inbred parents, whose inbreeding coefficient of 1
// is set rather than computed, but an inbreeding coefficient that comes out
// of the computation a rounding error short of 1, as for a self of a
// doubled haploid, leaves a variance of the order of 1e-16. Taken as a
// genuine variance, its reciprocal would put a coefficient of 1e16 into the
// inverse.
const double minVarD = 1e-10;

// A sparse row as (column, value) pairs
typedef std::vector<std::pair<arma::uword,double> > SpRow;

bool spRowLess(const std::pair<arma::uword,double>& a,
               const std::pair<arma::uword,double>& b){
  return a.first<b.first;
}

// Appends w times row to out, leaving repeated columns to tidyRow
void addScaledRow(SpRow& out, const SpRow& row, double w){
  for(SpRow::const_iterator it=row.begin(); it!=row.end(); ++it){
    out.push_back(std::make_pair(it->first, w*it->second));
  }
}

// Sorts a row by column, adds together entries in the same column and
// drops the ones that cancel
void tidyRow(SpRow& row){
  std::sort(row.begin(), row.end(), spRowLess);
  SpRow out;
  out.reserve(row.size());
  for(SpRow::const_iterator it=row.begin(); it!=row.end(); ++it){
    if(!out.empty() && (out.back().first==it->first)){
      out.back().second += it->second;
    }else{
      out.push_back(*it);
    }
  }
  SpRow kept;
  kept.reserve(out.size());
  for(SpRow::const_iterator it=out.begin(); it!=out.end(); ++it){
    if(it->second!=0.0){
      kept.push_back(*it);
    }
  }
  row.swap(kept);
}

// An individual's inbreeding coefficient from its parents' and the
// relationship between them, for ploidy v. Of the C(v,2) pairs of
// distinct alleles in the individual, C(v/2,2) lie within each parent's
// gamete and are identical by descent with that parent's inbreeding
// coefficient, and (v/2)^2 pair an allele from each gamete, identical with
// probability asd/v. For a diploid this is asd/2. An unknown parent is
// given F=0 and asd=0, as a gamete from the unrelated, non-inbred base.
double offspringF(int v, double Fs, double Fd, double asd){
  double h = 0.5*v;
  double within = 0.5*h*(h-1.0);
  double pairs = 0.5*v*(v-1.0);
  return (within*(Fs+Fd) + h*h*asd/v)/pairs;
}

// The variance of individual i's Mendelian sampling term, given the
// inbreeding coefficients of its parents. The variances below hold at
// every ploidy: writing the diagonal as 1+(v-1)F with F from offspringF,
// everything that depends on v cancels.
double mendelianVar(const Ped& ped, const arma::vec& F, arma::uword i){
  int s = ped.mother[i];
  int d = ped.father[i];
  if(ped.isDH[i]){
    // A doubled haploid is fully inbred, so its total variance is the
    // ploidy, 2. Its parent's breeding value accounts for 1+F of that, and
    // the rest is its own. An inbred founder has no parent and all of its
    // total variance is its own.
    if(s<0){
      return double(ped.ploidy);
    }
    return 1.0-F(s);
  }
  if((s>=0) && (d>=0)){
    return 0.5-0.25*(F(s)+F(d));
  }
  if(s>=0){
    return 0.75-0.25*F(s);
  }
  if(d>=0){
    return 0.75-0.25*F(d);
  }
  return 1.0;
}

// The column of the relationship matrix for individual p, over rows 0 to
// lim-1, by the two passes of Colleau's method; see calcPedASub. Every
// individual before lim must have its Mendelian sampling variance in D,
// and p must be before lim.
void colleauColumn(const Ped& ped, const arma::vec& D, arma::uword p,
                   arma::uword lim, arma::vec& x){
  x.zeros(lim);
  x(p) = 1.0;
  // Only p and its ancestors are nonzero on the way back
  for(arma::uword i=p+1; i-->0;){
    if(x(i)!=0.0){
      if(ped.mother[i]>=0){
        x(ped.mother[i]) += 0.5*x(i);
      }
      if(ped.father[i]>=0){
        x(ped.father[i]) += 0.5*x(i);
      }
    }
  }
  for(arma::uword i=0; i<lim; ++i){
    x(i) *= D(i);
    if(ped.mother[i]>=0){
      x(i) += 0.5*x(ped.mother[i]);
    }
    if(ped.father[i]>=0){
      x(i) += 0.5*x(ped.father[i]);
    }
  }
}

// Fills in F for the individuals in pend, all in one batch whose parents
// all come before the batch. Each needs half the relationship between its
// parents. The parents on whichever side has fewer distinct individuals
// each get one column of the relationship matrix, which gives the
// relationship to every partner at once.
void resolveBatch(const Ped& ped, std::vector<arma::uword>& pend,
                  const arma::vec& D, arma::vec& F, arma::vec& x){
  std::vector<int> mothers, fathers;
  for(std::vector<arma::uword>::const_iterator it=pend.begin();
      it!=pend.end(); ++it){
    mothers.push_back(ped.mother[*it]);
    fathers.push_back(ped.father[*it]);
  }
  std::sort(mothers.begin(), mothers.end());
  std::sort(fathers.begin(), fathers.end());
  arma::uword nMother = std::unique(mothers.begin(), mothers.end())-mothers.begin();
  arma::uword nFather = std::unique(fathers.begin(), fathers.end())-fathers.begin();
  const std::vector<int>& key = (nMother<=nFather) ? ped.mother : ped.father;
  const std::vector<int>& partner = (nMother<=nFather) ? ped.father : ped.mother;
  // Group the batch by the parent that gets the column
  std::sort(pend.begin(), pend.end(),
            [&key](arma::uword a, arma::uword c){ return key[a]<key[c]; });
  arma::uword g = 0;
  while(g<pend.size()){
    arma::uword p = key[pend[g]];
    arma::uword h = g;
    arma::uword lim = p+1;
    while((h<pend.size()) && (arma::uword(key[pend[h]])==p)){
      lim = std::max(lim, arma::uword(partner[pend[h]])+1);
      ++h;
    }
    colleauColumn(ped, D, p, lim, x);
    for(arma::uword k=g; k<h; ++k){
      F(pend[k]) = offspringF(ped.ploidy, F(key[pend[k]]),
                              F(partner[pend[k]]), x(partner[pend[k]]));
    }
    g = h;
  }
}

// Inbreeding coefficients already known, with NA for those that are not.
// An empty vector means that none are known.
arma::vec readKnownF(const Rcpp::NumericVector& knownF, arma::uword n){
  if((knownF.size()!=0) && (arma::uword(knownF.size())!=n)){
    Rcpp::stop("knownF must be empty or have one value per individual");
  }
  return Rcpp::as<arma::vec>(knownF);
}

// Returns a sparse matrix as 1-based triplets for R, which has no sparse
// matrix class without the Matrix package
Rcpp::List spTriplet(const arma::sp_mat& X){
  arma::uword nnz = X.n_nonzero;
  Rcpp::IntegerVector i(nnz), j(nnz);
  Rcpp::NumericVector x(nnz);
  arma::uword k = 0;
  for(arma::sp_mat::const_iterator it=X.begin(); it!=X.end(); ++it){
    i[k] = int(it.row())+1;
    j[k] = int(it.col())+1;
    x[k] = *it;
    ++k;
  }
  return Rcpp::List::create(Rcpp::Named("i")=i,
                            Rcpp::Named("j")=j,
                            Rcpp::Named("x")=x,
                            Rcpp::Named("nrow")=int(X.n_rows),
                            Rcpp::Named("ncol")=int(X.n_cols));
}

}

Ped readPed(const Rcpp::IntegerVector& mother,
            const Rcpp::IntegerVector& father,
            const Rcpp::IntegerVector& isDH,
            int ploidy){
  arma::uword n = mother.size();
  if((arma::uword(father.size())!=n) || (arma::uword(isDH.size())!=n)){
    Rcpp::stop("mother, father and isDH must have the same length");
  }
  if((ploidy<2) || ((ploidy%2)!=0)){
    Rcpp::stop("Pedigree relationships need an even ploidy level");
  }
  Ped ped;
  ped.n = n;
  ped.ploidy = ploidy;
  ped.mother.resize(n);
  ped.father.resize(n);
  ped.isDH.resize(n);
  for(arma::uword i=0; i<n; ++i){
    int m = mother[i];
    int f = father[i];
    int dh = isDH[i];
    // NA_INTEGER is the most negative int, so missing values fail the
    // first test. A 1-based parent of at most i is a 0-based row before i.
    if((m<0) || (f<0) || (m>int(i)) || (f>int(i))){
      Rcpp::stop("Every parent must come before its progeny in the pedigree");
    }
    if((dh!=0) && (dh!=1)){
      Rcpp::stop("isDH must be 0 or 1");
    }
    if((dh==1) && (m!=f)){
      Rcpp::stop("A doubled haploid must have the same mother and father");
    }
    if((dh==1) && (m>0) && (ploidy!=2)){
      Rcpp::stop("Doubled haploids are only supported for diploids");
    }
    ped.mother[i] = m-1;
    ped.father[i] = f-1;
    ped.isDH[i] = (dh==1);
  }
  return ped;
}

void pedInbreeding(const Ped& ped, const arma::vec& knownF, arma::vec& F,
                   arma::vec& D){
  arma::uword n = ped.n;
  bool useKnown = (knownF.n_elem>0);
  F.zeros(n);
  D.zeros(n);
  // The pedigree is taken in batches whose parents all come before the
  // batch, which is how AlphaSimR records the individuals made by one call.
  // Within a batch, an individual's inbreeding needs only the relationship
  // between its parents, and every relationship among earlier individuals
  // can be found by Colleau's passes once their variances are known. A
  // batch closes at the first individual with a parent inside it.
  std::vector<arma::uword> pend;
  arma::vec x;
  arma::uword b = 0;
  for(arma::uword i=0; i<=n; ++i){
    if((i==n) || (ped.mother[i]>=int(b)) || (ped.father[i]>=int(b))){
      if(i>b){
        if(!pend.empty()){
          resolveBatch(ped, pend, D, F, x);
          pend.clear();
        }
        for(arma::uword j=b; j<i; ++j){
          D(j) = mendelianVar(ped, F, j);
        }
      }
      b = i;
    }
    if(i==n){
      break;
    }
    int s = ped.mother[i];
    int d = ped.father[i];
    if(useKnown && !std::isnan(knownF(i))){
      F(i) = knownF(i);
    }else if(ped.isDH[i]){
      F(i) = 1.0;
    }else if((s<0) || (d<0)){
      // The unknown parent's gamete comes from the unrelated, non-inbred
      // base, so only the known parent's gamete can carry inbreeding,
      // which needs a ploidy above 2
      if((s<0) && (d<0)){
        F(i) = 0.0;
      }else{
        F(i) = offspringF(ped.ploidy, (s>=0) ? F(s) : F(d), 0.0, 0.0);
      }
    }else if(s==d){
      // A self needs only its parent's own relationship, 1+(v-1)F
      F(i) = offspringF(ped.ploidy, F(s), F(s),
                        1.0+(ped.ploidy-1.0)*F(s));
    }else{
      pend.push_back(i);
    }
  }
}

PedAinv buildPedAinv(const Ped& ped, const arma::vec& D){
  arma::uword n = ped.n;
  // w[i] expresses individual i's breeding value in terms of the basis
  std::vector<SpRow> w(n);
  std::vector<arma::uword> basis;
  basis.reserve(n);
  std::vector<arma::uword> ti, tj;
  std::vector<double> tx;
  SpRow c;
  for(arma::uword i=0; i<n; ++i){
    int s = ped.mother[i];
    int d = ped.father[i];
    // The parents' share of i's breeding value
    c.clear();
    if(s>=0){
      addScaledRow(c, w[s], 0.5);
    }
    if(d>=0){
      addScaledRow(c, w[d], 0.5);
    }
    if(D(i)>minVarD){
      arma::uword k = basis.size();
      basis.push_back(i);
      w[i].assign(1, std::make_pair(k, 1.0));
      // i's Mendelian sampling term is its breeding value less its parents'
      // share, with variance D(i). Its contribution to the inverse is the
      // outer product of the coefficients of that difference, divided by
      // D(i). With known parents in the basis, this is Henderson's rule of
      // 1, -1/2 and 1/4 over 1+F terms. A doubled haploid's parent gets
      // -1/2 twice, and 1/4 four times, which sum to the -1 and 1 of a
      // single parent with the whole share.
      for(SpRow::iterator it=c.begin(); it!=c.end(); ++it){
        it->second = -it->second;
      }
      c.push_back(std::make_pair(k, 1.0));
      tidyRow(c);
      double b = 1.0/D(i);
      for(SpRow::const_iterator it1=c.begin(); it1!=c.end(); ++it1){
        for(SpRow::const_iterator it2=c.begin(); it2!=c.end(); ++it2){
          ti.push_back(it1->first);
          tj.push_back(it2->first);
          tx.push_back(b*it1->second*it2->second);
        }
      }
    }else{
      // No Mendelian sampling, so i is exactly its parents' share
      tidyRow(c);
      w[i] = c;
    }
  }

  arma::uword nBasis = basis.size();
  PedAinv output;
  output.basis = arma::conv_to<arma::uvec>::from(basis);

  arma::umat loc(2, tx.size());
  arma::vec val(tx.size());
  for(arma::uword k=0; k<tx.size(); ++k){
    loc(0,k) = ti[k];
    loc(1,k) = tj[k];
    val(k) = tx[k];
  }
  output.Ainv = arma::sp_mat(true, loc, val, nBasis, nBasis);

  arma::uword nW = 0;
  for(arma::uword i=0; i<n; ++i){
    nW += w[i].size();
  }
  arma::umat wLoc(2, nW);
  arma::vec wVal(nW);
  arma::uword k = 0;
  for(arma::uword i=0; i<n; ++i){
    for(SpRow::const_iterator it=w[i].begin(); it!=w[i].end(); ++it){
      wLoc(0,k) = i;
      wLoc(1,k) = it->first;
      wVal(k) = it->second;
      ++k;
    }
  }
  output.W = arma::sp_mat(true, wLoc, wVal, n, nBasis);

  return output;
}

// Inbreeding coefficients of a pedigree, see readPed for its form. knownF
// holds values already computed, with NA for the rest.
// [[Rcpp::export]]
arma::vec calcPedInbreeding(const Rcpp::IntegerVector& mother,
                            const Rcpp::IntegerVector& father,
                            const Rcpp::IntegerVector& isDH,
                            const Rcpp::NumericVector& knownF=Rcpp::NumericVector::create(),
                            int ploidy=2){
  Ped ped = readPed(mother, father, isDH, ploidy);
  arma::vec F, D;
  pedInbreeding(ped, readKnownF(knownF, ped.n), F, D);
  return F;
}

// The inverse of the relationship matrix as triplets, with the matrix W
// that maps the basis individuals onto the whole pedigree. See PedAinv in
// pedigree.h.
// [[Rcpp::export]]
Rcpp::List calcPedAinvTriplet(const Rcpp::IntegerVector& mother,
                              const Rcpp::IntegerVector& father,
                              const Rcpp::IntegerVector& isDH,
                              const Rcpp::NumericVector& knownF=Rcpp::NumericVector::create(),
                              int ploidy=2){
  Ped ped = readPed(mother, father, isDH, ploidy);
  arma::vec F, D;
  pedInbreeding(ped, readKnownF(knownF, ped.n), F, D);
  PedAinv ans = buildPedAinv(ped, D);
  Rcpp::IntegerVector basis(ans.basis.n_elem);
  for(arma::uword k=0; k<ans.basis.n_elem; ++k){
    basis[k] = int(ans.basis(k))+1;
  }
  return Rcpp::List::create(Rcpp::Named("Ainv")=spTriplet(ans.Ainv),
                            Rcpp::Named("W")=spTriplet(ans.W),
                            Rcpp::Named("basis")=basis,
                            Rcpp::Named("F")=F,
                            Rcpp::Named("D")=D);
}

// The relationship matrix by the tabular method. Each column is the
// average of the parents' columns, so it is filled one individual at a
// time in pedigree order, writing the row as well as the column so that
// every parent's column is complete when it is read. It needs the whole
// matrix in memory, so it is meant for pedigrees of modest size.
// [[Rcpp::export]]
arma::mat calcPedA(const Rcpp::IntegerVector& mother,
                   const Rcpp::IntegerVector& father,
                   const Rcpp::IntegerVector& isDH,
                   int ploidy=2){
  Ped ped = readPed(mother, father, isDH, ploidy);
  arma::uword n = ped.n;
  arma::mat A(n, n, arma::fill::zeros);
  for(arma::uword i=0; i<n; ++i){
    int s = ped.mother[i];
    int d = ped.father[i];
    for(arma::uword j=0; j<i; ++j){
      double a = 0.0;
      if(s>=0){
        a += 0.5*A(j,s);
      }
      if(d>=0){
        a += 0.5*A(j,d);
      }
      A(j,i) = a;
      A(i,j) = a;
    }
    double v = ped.ploidy;
    if(ped.isDH[i]){
      // An inbred founder, or a diploid's doubled haploid
      A(i,i) = v;
    }else if((s>=0) || (d>=0)){
      double Fs = (s>=0) ? (A(s,s)-1.0)/(v-1.0) : 0.0;
      double Fd = (d>=0) ? (A(d,d)-1.0)/(v-1.0) : 0.0;
      double asd = ((s>=0) && (d>=0)) ? A(s,d) : 0.0;
      A(i,i) = 1.0+(v-1.0)*offspringF(ped.ploidy, Fs, Fd, asd);
    }else{
      A(i,i) = 1.0;
    }
  }
  return A;
}

// The block of the relationship matrix among the target individuals, by
// the indirect method of Colleau (2002, doi:10.1186/1297-9686-34-4-409).
// A = T*diag(D)*T', and T's inverse is the identity less half of each
// parent, so the column of A for one individual is found with two passes
// over the pedigree: one from the latest individual back to the founders,
// which passes half of each individual's value to each of its parents and
// gives T'*x, and one forward from the founders, after scaling by D, which
// adds half of each parent's value to its progeny and gives T*(D*(T'*x)).
// A doubled haploid's parent receives both halves, as everywhere else in
// this file. D is only multiplied, never inverted, so individuals with no
// Mendelian sampling need no special handling here. Memory use is that of
// the block plus one vector over the pedigree, and the time that of one
// pass over the pedigree per target, once the inbreeding coefficients are
// known. knownF holds those already computed, with NA for the rest.
// [[Rcpp::export]]
arma::mat calcPedASub(const Rcpp::IntegerVector& mother,
                      const Rcpp::IntegerVector& father,
                      const Rcpp::IntegerVector& isDH,
                      const Rcpp::IntegerVector& target,
                      const Rcpp::NumericVector& knownF=Rcpp::NumericVector::create(),
                      int ploidy=2){
  Ped ped = readPed(mother, father, isDH, ploidy);
  arma::uword n = ped.n;
  arma::uword k = target.size();
  std::vector<arma::uword> take(k);
  for(arma::uword c=0; c<k; ++c){
    int t = target[c];
    // NA_INTEGER is the most negative int, so it fails this test too
    if((t<1) || (t>int(n))){
      Rcpp::stop("Every target must be a row of the pedigree");
    }
    take[c] = arma::uword(t-1);
  }
  arma::vec F, D;
  pedInbreeding(ped, readKnownF(knownF, ped.n), F, D);
  arma::mat output(k, k);
  arma::vec x(n);
  for(arma::uword c=0; c<k; ++c){
    x.zeros();
    x(take[c]) = 1.0;
    // Only the target and its ancestors are nonzero on the way back
    for(arma::uword i=n; i-->0;){
      if(x(i)!=0.0){
        if(ped.mother[i]>=0){
          x(ped.mother[i]) += 0.5*x(i);
        }
        if(ped.father[i]>=0){
          x(ped.father[i]) += 0.5*x(i);
        }
      }
    }
    x %= D;
    for(arma::uword i=0; i<n; ++i){
      if(ped.mother[i]>=0){
        x(i) += 0.5*x(ped.mother[i]);
      }
      if(ped.father[i]>=0){
        x(i) += 0.5*x(ped.father[i]);
      }
    }
    for(arma::uword r=0; r<k; ++r){
      output(r,c) = x(take[r]);
    }
  }
  // Each column is computed separately, so rounding can leave the two
  // halves a few units in the last place apart
  output = 0.5*(output+output.t());
  return output;
}

// The parts of the mixed model equations for pedigree BLUP of one trait,
// as triplets for Matrix's sparse Cholesky in R. The model is
// y = X*b + Z*W*u + e, where X has one column per level of the fixed
// effects, Z maps each record to its row of the pedigree, and W maps the
// pedigree onto the basis of PedAinv, so that u ~ N(0, Ainv^-1*Vu) over the
// basis. The coefficient matrix is C0 + lambda*A0 with lambda = Ve/Vu,
// where C0 holds X'X, X'ZW and W'Z'ZW, and A0 the inverse relationship
// matrix, each as an upper triangle with the fixed effects first and
// repeated entries to be summed. rec and fix are 1-based and given for
// each record, and the levels of fix must run from 1 with none missing.
// [[Rcpp::export]]
Rcpp::List calcPedMME(const Rcpp::IntegerVector& mother,
                      const Rcpp::IntegerVector& father,
                      const Rcpp::IntegerVector& isDH,
                      const Rcpp::IntegerVector& rec,
                      const Rcpp::IntegerVector& fix,
                      const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& knownF=Rcpp::NumericVector::create(),
                      int ploidy=2){
  Ped ped = readPed(mother, father, isDH, ploidy);
  arma::uword N = y.size();
  if((arma::uword(rec.size())!=N) || (arma::uword(fix.size())!=N)){
    Rcpp::stop("rec, fix and y must have the same length");
  }
  int nFix = 0;
  for(arma::uword r=0; r<N; ++r){
    if((rec[r]<1) || (rec[r]>int(ped.n))){
      Rcpp::stop("Every record must belong to a row of the pedigree");
    }
    if(fix[r]<1){
      Rcpp::stop("Fixed effect levels must start from 1");
    }
    if(!std::isfinite(y[r])){
      Rcpp::stop("Records must be finite");
    }
    nFix = std::max(nFix, int(fix[r]));
  }
  arma::vec F, D;
  pedInbreeding(ped, readKnownF(knownF, ped.n), F, D);
  PedAinv ans = buildPedAinv(ped, D);
  arma::uword nB = ans.basis.n_elem;

  // The rows of W, which Armadillo stores by column
  std::vector<SpRow> w(ped.n);
  for(arma::sp_mat::const_iterator it=ans.W.begin(); it!=ans.W.end(); ++it){
    w[it.row()].push_back(std::make_pair(arma::uword(it.col()), double(*it)));
  }

  std::vector<int> ci, cj;
  std::vector<double> cx;
  std::vector<double> nRec(nFix, 0.0);
  arma::vec rhs(nFix+nB, arma::fill::zeros);
  double yy = 0.0;
  for(arma::uword r=0; r<N; ++r){
    arma::uword f = fix[r]-1;
    const SpRow& row = w[rec[r]-1];
    double yr = y[r];
    nRec[f] += 1.0;
    rhs(f) += yr;
    yy += yr*yr;
    // A record of a basis individual touches one random effect, and one of
    // a clone or an F1 of inbred parents touches its basis ancestors
    for(SpRow::const_iterator e1=row.begin(); e1!=row.end(); ++e1){
      ci.push_back(int(f)+1);
      cj.push_back(nFix+int(e1->first)+1);
      cx.push_back(e1->second);
      rhs(nFix+e1->first) += e1->second*yr;
      for(SpRow::const_iterator e2=row.begin(); e2!=row.end(); ++e2){
        if(e1->first<=e2->first){
          ci.push_back(nFix+int(e1->first)+1);
          cj.push_back(nFix+int(e2->first)+1);
          cx.push_back(e1->second*e2->second);
        }
      }
    }
  }
  for(int f=0; f<nFix; ++f){
    ci.push_back(f+1);
    cj.push_back(f+1);
    cx.push_back(nRec[f]);
  }

  std::vector<int> ai, aj;
  std::vector<double> ax;
  for(arma::sp_mat::const_iterator it=ans.Ainv.begin(); it!=ans.Ainv.end(); ++it){
    if(it.row()<=it.col()){
      ai.push_back(nFix+int(it.row())+1);
      aj.push_back(nFix+int(it.col())+1);
      ax.push_back(*it);
    }
  }

  // log|A| over the basis, the constant of the REML likelihood that comes
  // from the relationships, is the sum of the Mendelian sampling variances'
  // logarithms
  double logDetA = 0.0;
  for(arma::uword k=0; k<nB; ++k){
    logDetA += std::log(D(ans.basis(k)));
  }

  Rcpp::IntegerVector basis(nB);
  for(arma::uword k=0; k<nB; ++k){
    basis[k] = int(ans.basis(k))+1;
  }
  return Rcpp::List::create(
    Rcpp::Named("C")=Rcpp::List::create(Rcpp::Named("i")=ci,
                                        Rcpp::Named("j")=cj,
                                        Rcpp::Named("x")=cx),
    Rcpp::Named("A")=Rcpp::List::create(Rcpp::Named("i")=ai,
                                        Rcpp::Named("j")=aj,
                                        Rcpp::Named("x")=ax),
    Rcpp::Named("rhs")=rhs,
    Rcpp::Named("yy")=yy,
    Rcpp::Named("nFix")=nFix,
    Rcpp::Named("nBasis")=int(nB),
    Rcpp::Named("logDetA")=logDetA,
    Rcpp::Named("W")=spTriplet(ans.W),
    Rcpp::Named("basis")=basis);
}
