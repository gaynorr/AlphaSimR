#include "alphasimr.h"
#include "pedigree.h"
#include <algorithm>
#include <queue>
#include <utility>

// Numerator relationships from a pedigree: inbreeding coefficients by the
// method of Meuwissen and Luo (1992, doi:10.1186/1297-9686-24-4-305), the
// inverse of the relationship matrix by the rules of Henderson (1976,
// doi:10.2307/2529339) and Quaas (1976, doi:10.2307/2529279), and the
// relationship matrix itself by Henderson's tabular method.
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

// The variance of individual i's Mendelian sampling term, given the
// inbreeding coefficients of its parents
double mendelianVar(const Ped& ped, const arma::vec& F, arma::uword i){
  int s = ped.mother[i];
  int d = ped.father[i];
  if(ped.isDH[i]){
    // A doubled haploid is fully inbred, so its total variance is 2. Its
    // parent's breeding value accounts for 1+F of that, and the rest is its
    // own. An inbred founder has no parent and all of the 2 is its own.
    if(s<0){
      return 2.0;
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
            const Rcpp::IntegerVector& isDH){
  arma::uword n = mother.size();
  if((arma::uword(father.size())!=n) || (arma::uword(isDH.size())!=n)){
    Rcpp::stop("mother, father and isDH must have the same length");
  }
  Ped ped;
  ped.n = n;
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
    ped.mother[i] = m-1;
    ped.father[i] = f-1;
    ped.isDH[i] = (dh==1);
  }
  return ped;
}

void pedInbreeding(const Ped& ped, arma::vec& F, arma::vec& D){
  arma::uword n = ped.n;
  F.zeros(n);
  D.zeros(n);
  // Row i of T holds the share of each ancestor's Mendelian sampling term
  // in individual i, and A(i,i) is the sum of those shares squared times
  // the variances. The row is filled from i back toward the founders,
  // always taking the latest ancestor next, so that every contribution to
  // an ancestor's share has arrived before it is used. Parents come before
  // their progeny, so an ancestor that has been taken is never reached
  // again.
  std::vector<double> T(n, 0.0);
  std::priority_queue<arma::uword> next;
  for(arma::uword i=0; i<n; ++i){
    int s = ped.mother[i];
    int d = ped.father[i];
    D(i) = mendelianVar(ped, F, i);
    if(ped.isDH[i]){
      F(i) = 1.0;
    }else if((s<0) || (d<0)){
      // With at most one known parent, no two alleles can be traced back
      // to a common ancestor
      F(i) = 0.0;
    }else if((i>0) && (s==ped.mother[i-1]) && (d==ped.father[i-1]) &&
             !ped.isDH[i-1]){
      // Full sibs have the same inbreeding coefficient, and AlphaSimR
      // records the progeny of a cross together
      F(i) = F(i-1);
    }else{
      double aii = 0.0;
      T[i] = 1.0;
      next.push(i);
      while(!next.empty()){
        arma::uword j = next.top();
        next.pop();
        double t = T[j];
        T[j] = 0.0;
        aii += t*t*D(j);
        int sj = ped.mother[j];
        int dj = ped.father[j];
        // Shares are always positive, so a share of zero means that the
        // ancestor has not been queued yet
        if(sj>=0){
          if(T[sj]==0.0){
            next.push(sj);
          }
          T[sj] += 0.5*t;
        }
        if(dj>=0){
          if(T[dj]==0.0){
            next.push(dj);
          }
          T[dj] += 0.5*t;
        }
      }
      F(i) = aii-1.0;
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

// Inbreeding coefficients of a pedigree, see readPed for its form
// [[Rcpp::export]]
arma::vec calcPedInbreeding(const Rcpp::IntegerVector& mother,
                            const Rcpp::IntegerVector& father,
                            const Rcpp::IntegerVector& isDH){
  Ped ped = readPed(mother, father, isDH);
  arma::vec F, D;
  pedInbreeding(ped, F, D);
  return F;
}

// The inverse of the relationship matrix as triplets, with the matrix W
// that maps the basis individuals onto the whole pedigree. See PedAinv in
// pedigree.h.
// [[Rcpp::export]]
Rcpp::List calcPedAinvTriplet(const Rcpp::IntegerVector& mother,
                              const Rcpp::IntegerVector& father,
                              const Rcpp::IntegerVector& isDH){
  Ped ped = readPed(mother, father, isDH);
  arma::vec F, D;
  pedInbreeding(ped, F, D);
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
                   const Rcpp::IntegerVector& isDH){
  Ped ped = readPed(mother, father, isDH);
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
    if(ped.isDH[i]){
      A(i,i) = 2.0;
    }else if((s>=0) && (d>=0)){
      A(i,i) = 1.0+0.5*A(s,d);
    }else{
      A(i,i) = 1.0;
    }
  }
  return A;
}
