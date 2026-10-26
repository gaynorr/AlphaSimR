#include "alphasimr.h"
#include "packedGv.h"

// Eigen decomposition of a small symmetric matrix by cyclic Jacobi
// rotations. It is used inside a parallel loop on matrices of at most four
// rows, where calling LAPACK is unnecessary and its thread safety depends on
// the BLAS the user has installed. The result is deterministic.
static void jacobiEigen(arma::mat S, arma::vec& lambda, arma::mat& V){
  const arma::uword n = S.n_rows;
  V.eye(n, n);
  for(int sweep=0; sweep<100; ++sweep){
    double off = 0.0, total = 0.0;
    for(arma::uword p=0; p<n; ++p){
      for(arma::uword q=0; q<n; ++q){
        total += S(p,q)*S(p,q);
        if(p!=q) off += S(p,q)*S(p,q);
      }
    }
    if(off<=1e-30*total){
      break;
    }
    for(arma::uword p=0; p+1<n; ++p){
      for(arma::uword q=p+1; q<n; ++q){
        if(S(p,q)==0.0) continue;
        double theta = (S(q,q)-S(p,p))/(2.0*S(p,q));
        double t = (theta>=0.0 ? 1.0 : -1.0)/
          (std::abs(theta)+std::sqrt(theta*theta+1.0));
        double c = 1.0/std::sqrt(t*t+1.0);
        double s = t*c;
        for(arma::uword k=0; k<n; ++k){
          double skp = S(k,p), skq = S(k,q);
          S(k,p) = c*skp-s*skq;
          S(k,q) = s*skp+c*skq;
        }
        for(arma::uword k=0; k<n; ++k){
          double spk = S(p,k), sqk = S(q,k);
          S(p,k) = c*spk-s*sqk;
          S(q,k) = s*spk+c*sqk;
        }
        for(arma::uword k=0; k<n; ++k){
          double vkp = V(k,p), vkq = V(k,q);
          V(k,p) = c*vkp-s*vkq;
          V(k,q) = s*vkp+c*vkq;
        }
      }
    }
  }
  lambda = S.diag();
}

// Weighted least squares of y on the first nFit columns of Z, with weights
// w, giving the minimum norm solution. The columns are first scaled to unit
// weighted sum of squares, and a column whose sum of squares is negligible,
// such as the dosage of a fixed locus, gets a coefficient of zero. The
// scaled normal equations are solved through their eigen decomposition,
// leaving out directions the data do not identify. The fitted values are
// therefore the projection of y onto the span of the columns even when the
// fit is singular, as it is when the two loci of a pair are in perfect
// linkage disequilibrium. The coefficients are then not identified, and the
// minimum norm solution splits them evenly between columns that carry the
// same information, so the result does not depend on which locus of a pair
// comes first. The systems are tiny and are solved inside a parallel loop,
// which is why no library routine is called.
static arma::vec wlsFit(const arma::mat& Z, const arma::vec& w,
                        const arma::vec& y, arma::uword nFit){
  const double tolAbs = 1e-14, tolRel = 1e-10;
  arma::vec coef(nFit, arma::fill::zeros);
  std::vector<arma::uword> keep;
  arma::vec scale(nFit, arma::fill::zeros);
  for(arma::uword c=0; c<nFit; ++c){
    double ss = arma::accu(w%Z.col(c)%Z.col(c));
    if(ss>tolAbs){
      keep.push_back(c);
      scale(c) = std::sqrt(ss);
    }
  }
  arma::uword k = keep.size();
  if(k==0){
    return coef;
  }
  arma::mat S(k, k);
  arma::vec r(k);
  for(arma::uword p=0; p<k; ++p){
    arma::vec zp = Z.col(keep[p])/scale(keep[p]);
    r(p) = arma::accu(w%zp%y);
    for(arma::uword q=0; q<=p; ++q){
      arma::vec zq = Z.col(keep[q])/scale(keep[q]);
      S(p,q) = arma::accu(w%zp%zq);
      S(q,p) = S(p,q);
    }
  }
  arma::vec lambda;
  arma::mat V;
  jacobiEigen(S, lambda, V);
  double lambdaMax = lambda.max();
  arma::vec b(k, arma::fill::zeros);
  for(arma::uword m=0; m<k; ++m){
    if(lambda(m)>tolRel*lambdaMax){
      double proj = 0.0;
      for(arma::uword p=0; p<k; ++p){
        proj += V(p,m)*r(p);
      }
      proj /= lambda(m);
      for(arma::uword p=0; p<k; ++p){
        b(p) += V(p,m)*proj;
      }
    }
  }
  for(arma::uword p=0; p<k; ++p){
    coef(keep[p]) = b(p)/scale(keep[p]);
  }
  return coef;
}

// Calculates genetic parameters for traits with epistasis
//
// Genotypes are not decoded into a matrix. One pass over tiles of
// individuals counts the joint genotypes of every pair of loci and adds up
// each individual's epistatic genetic value, gv_aa. The statistics for each
// pair follow from its counts. Every locus belongs to exactly one pair, and
// once a pair's frequencies are known a locus's breeding value and
// dominance deviation are functions of its own dosage with the same form as
// a trait's additive and dominance effects plus a constant, so gv_a, gv_d,
// bv and dd are sets of effects on the trait's loci and are worked out by
// gvGroupEffects. The epistatic deviation is what is left of an
// individual's genetic value: gv_a+gv_d+gv_aa-bv-dd less the sum of each
// pair's mean genetic value.
//
// A locus's own additive and dominance effects are regressed on its own
// dosage, as for a trait without epistasis, so linkage disequilibrium
// between loci appears in the covariances between them. The epistatic value
// of a pair is a function of both dosages, so it is fitted to the two
// dosages jointly, under the pair's observed joint genotype frequencies, by
// weighted least squares: the fit on the two dosages gives its part of the
// breeding values, and the further fit on the two dominance codings gives
// its part of the dominance deviations. Under linkage disequilibrium this is
// what makes epistatic variance additive or dominance variance, and under
// linkage equilibrium it reduces to averaging over the partner locus. The
// genic variances, and the tables for genic covariances, are calculated
// under linkage equilibrium and do not use the joint fit.
//
// With locusTables, the genotype frequencies and statistical values at each
// locus, and the epistatic deviations of each pair, are also returned, so
// that genParam can form genic covariances between traits.
Rcpp::List calcGenParamE(const Rcpp::S4& trait,
                         const Rcpp::S4& pop,
                         int nThreads,
                         bool locusTables=false){
  //Information from pop
  bool hasD = trait.hasSlot("domEff");
  PackedGeno g(pop);
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  double dP = double(ploidy);
  arma::uword nGeno = ploidy+1;
  //Information from trait
  const arma::Col<int>& lociPerChr = trait.slot("lociPerChr");
  arma::uvec lociLoc = trait.slot("lociLoc");
  arma::vec a = trait.slot("addEff");
  arma::mat E;
  E = Rcpp::as<arma::mat>(trait.slot("epiEff"));
  E.col(0) -= 1; //R to C++
  E.col(1) -= 1; //R to C++
  arma::vec d;
  if(hasD){
    d = Rcpp::as<arma::vec>(trait.slot("domEff"));
  }
  double intercept = trait.slot("intercept");
  arma::uword nLoci = a.n_elem;
  arma::uword nPair = E.n_rows;
  LociPos pos = findLoci(g, lociPerChr, lociLoc);
  if(pos.base.size()!=nLoci){
    Rcpp::stop("number of loci does not match the number of effects");
  }
  arma::uvec loc1(nPair), loc2(nPair);
  for(arma::uword i=0; i<nPair; ++i){
    if(E(i,0)<0 || E(i,1)<0 || E(i,0)>=nLoci || E(i,1)>=nLoci){
      Rcpp::stop("epiEff refers to a locus the trait does not have");
    }
    loc1(i) = arma::uword(E(i,0));
    loc2(i) = arma::uword(E(i,1));
  }
  arma::vec x(ploidy+1); // Genotype dosage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);

  // Each pair's epistatic genetic value over its two dosages
  arma::mat aaTab(nGeno*nGeno, nPair);
  for(arma::uword i=0; i<nPair; ++i){
    arma::mat aaEff = xa*xa.t()*E(i,2);
    aaTab.col(i) = arma::vectorise(aaEff);
  }

  // Joint genotype counts and gv_aa. Individuals are split into blocks that
  // own their part of gv_aa, and each individual sums its pairs in order,
  // so gv_aa does not depend on the thread count. Counts are integers kept
  // by each thread and added afterward, which is exact in any order.
  arma::vec gvAA(nInd, arma::fill::zeros);
  arma::uword nIndBlocks = countBlocks(nInd);
  int nThreadsInd = nThreads;
  if(nIndBlocks < static_cast<arma::uword>(nThreadsInd)){
    nThreadsInd = static_cast<int>(nIndBlocks);
  }
  if(nThreadsInd<1){
    nThreadsInd = 1;
  }
  arma::uword nCell = nGeno*nGeno;
  std::vector<std::vector<uint32_t> > threadCounts(nThreadsInd,
                                                   std::vector<uint32_t>(nCell*nPair, 0));
  DosageTiles tiles(g, pos);
  const arma::uword kTile = 32;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreadsInd)
#endif
  for(arma::uword block=0; block<nIndBlocks; ++block){
#ifdef _OPENMP
    uint32_t* cnt = threadCounts[omp_get_thread_num()].data();
#else
    uint32_t* cnt = threadCounts[0].data();
#endif
    arma::uword indStart = blockStart(nInd, nIndBlocks, block);
    arma::uword indEnd = blockStart(nInd, nIndBlocks, block+1);
    std::vector<unsigned char> buf(nLoci*kTile);
    for(arma::uword t0=indStart; t0<indEnd; t0+=kTile){
      arma::uword nj = std::min(t0+kTile, indEnd)-t0;
      tiles.decode(t0, nj, kTile, buf.data());
      double* aaOut = gvAA.memptr()+t0;
      for(arma::uword i=0; i<nPair; ++i){
        const unsigned char* u = &buf[loc1(i)*kTile];
        const unsigned char* v = &buf[loc2(i)*kTile];
        const double* tab = aaTab.colptr(i);
        uint32_t* c = cnt+i*nCell;
        for(arma::uword j=0; j<nj; ++j){
          arma::uword uv = u[j]+v[j]*nGeno;
          aaOut[j] += tab[uv];
          ++c[uv];
        }
      }
    }
  }
  arma::Mat<arma::uword> pairCounts(nCell, nPair, arma::fill::zeros);
  for(int t=0; t<nThreadsInd; ++t){
    const uint32_t* c = threadCounts[t].data();
    arma::uword* pc = pairCounts.memptr();
    for(arma::uword k=0; k<nCell*nPair; ++k){
      pc[k] += c[k];
    }
  }
  threadCounts.clear();

  // The statistics for each pair. Accumulators are indexed by work block,
  // not by thread, so that their number and the order they are summed in
  // do not depend on how many threads are available
  arma::uword nBlocks = countBlocks(nPair);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::vec genicA(nBlocks,arma::fill::zeros); // No LD
  arma::vec genicA2(nBlocks,arma::fill::zeros); // No LD and HWE
  arma::vec genicD(nBlocks,arma::fill::zeros); // No LD
  arma::vec genicD2(nBlocks,arma::fill::zeros); // No LD and HWE
  arma::vec genicAA(nBlocks,arma::fill::zeros); // No LD
  arma::vec genicAA2(nBlocks,arma::fill::zeros); // No LD and HWE
  arma::vec mu(nBlocks,arma::fill::zeros); // Observed mean
  arma::vec eMu(nBlocks,arma::fill::zeros); // Expected mean with HWE
  arma::vec leMu(nBlocks,arma::fill::zeros); // Expected mean with no LD
  arma::vec alpha(a.n_elem);
  arma::vec alphaHW(a.n_elem);
  // For each locus: its mean dosage, and the additive and dominance effects
  // and constant that give its dominance deviation. For each pair: its mean
  // genetic value.
  arma::vec locGenoMu(nLoci), locDdAdd(nLoci), locDdDom(nLoci),
    locDdConst(nLoci), pairGvMu(nPair);
  // Tables for genic covariances between traits: one column per locus, or
  // per pair for the epistatic deviations. Each column is written by the
  // one pair that owns it, so the writes need no synchronization, and none
  // of these values feeds the sums above.
  arma::mat tabFreq, tabFreqE, tabBv, tabBvE, tabDd, tabDdE, tabAA, tabAAE;
  if(locusTables){
    tabFreq.set_size(nGeno, nLoci);
    tabFreqE.set_size(nGeno, nLoci);
    tabBv.set_size(nGeno, nLoci);
    tabBvE.set_size(nGeno, nLoci);
    if(hasD){
      tabDd.set_size(nGeno, nLoci);
      tabDdE.set_size(nGeno, nLoci);
    }
    tabAA.set_size(nCell, nPair);
    tabAAE.set_size(nCell, nPair);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword tid=0; tid<nBlocks; ++tid){
    arma::uword itemStart = blockStart(E.n_rows, nBlocks, tid);
    arma::uword itemEnd = blockStart(E.n_rows, nBlocks, tid+1);
    for(arma::uword i=itemStart; i<itemEnd; ++i){
      double gvMu1, gvMu2, gvEMu1, gvEMu2, 
      genoMu1, genoMu2, p1, p2, q1, q2, dK, 
      gvMu, gvEMu, gvNoLDMu;
    
    
      //Observed frequencies
      arma::mat freq(ploidy+1,ploidy+1);
      for(arma::uword v=0; v<(ploidy+1); ++v){
        for(arma::uword u=0; u<(ploidy+1); ++u){
          freq(u,v) = double(pairCounts(u+v*(ploidy+1),i));
        }
      }
      freq = freq/accu(freq);
      arma::vec freq1 = sum(freq,1);
      arma::vec freq2 = sum(freq,0).t();
    
      genoMu1 = accu(freq1%x);
      p1 = genoMu1/dP;
      q1 = 1-p1;
    
      genoMu2 = accu(freq2%x);
      p2 = genoMu2/dP;
      q2 = 1-p2;
    
      // Expected frequencies
      arma::vec freqE1(ploidy+1), freqE2(ploidy+1);
      for(arma::uword k=0; k<(ploidy+1); ++k){
        dK = double(k);
        freqE1(k) = choose(dP,dK)*std::pow(p1,dK)*std::pow(q1,dP-dK);
        freqE2(k) = choose(dP,dK)*std::pow(p2,dK)*std::pow(q2,dP-dK);
      }
    
      // Frequencies with no LD
      arma::mat freqNoLD(ploidy+1,ploidy+1);
      arma::mat freqNoLDE(ploidy+1,ploidy+1);
      for(arma::uword j=0; j<(ploidy+1); ++j){
        for(arma::uword k=0; k<(ploidy+1); ++k){
          freqNoLDE(j,k) = freqE1(j)*freqE2(k);
          freqNoLD(j,k) = freq1(j)*freq2(k);
        }
      }
    
      //Marginal values (individual loci)
      //Additive effects
      arma::vec aEff1 = xa*a(E(i,0));
      arma::vec aEff2 = xa*a(E(i,1));
      //Additive-by-additive effects
      arma::mat aaEff = xa*xa.t()*E(i,2);
      //Dominance effects
      arma::vec dEff1, dEff2;
      //Genetic value
      arma::vec gv1, gv2, gvE1, gvE2;
      if(hasD){
        dEff1 = xd*d(E(i,0));
        dEff2 = xd*d(E(i,1));
        gv1 = aEff1+dEff1;
        gv2 = aEff2+dEff2;
        gvE1 = gv1;
        gvE2 = gv2;
        for(arma::uword j=0; j<(ploidy+1); ++j){
          gv1(j) += accu(freq2%(aEff2+dEff2+E(i,2)*xa(j)*xa));
          gv2(j) += accu(freq1%(aEff1+dEff1+E(i,2)*xa(j)*xa));
          gvE1(j) += accu(freqE2%(aEff2+dEff2+E(i,2)*xa(j)*xa));
          gvE2(j) += accu(freqE1%(aEff1+dEff1+E(i,2)*xa(j)*xa));
        }
      }else{
        gv1 = aEff1;
        gv2 = aEff2;
        gvE1 = gv1;
        gvE2 = gv2;
        for(arma::uword j=0; j<(ploidy+1); ++j){
          gv1(j) += accu(freq2%(aEff2+E(i,2)*xa(j)*xa));
          gv2(j) += accu(freq1%(aEff1+E(i,2)*xa(j)*xa));
          gvE1(j) += accu(freqE2%(aEff2+E(i,2)*xa(j)*xa));
          gvE2(j) += accu(freqE1%(aEff1+E(i,2)*xa(j)*xa));
        }
      }
    
      gvMu1 = accu(freq1%gv1);
      gvMu2 = accu(freq2%gv2);
      gvEMu1 = accu(freqE1%gvE1);
      gvEMu2 = accu(freqE2%gvE2);
    
      // Average effects with the partner averaged over its marginal
      // frequencies, as under linkage equilibrium. These give the genic
      // variances and tables; the observed average effects are worked out
      // below from the joint frequencies.
      double alphaLe1 = accu(freq1%(gv1-gvMu1)%(x-genoMu1))/
        accu(freq1%(x-genoMu1)%(x-genoMu1));
      alphaHW(E(i,0)) = accu(freqE1%(gvE1-gvEMu1)%(x-genoMu1))/
        accu(freqE1%(x-genoMu1)%(x-genoMu1));
      double alphaLe2 = accu(freq2%(gv2-gvMu2)%(x-genoMu2))/
        accu(freq2%(x-genoMu2)%(x-genoMu2));
      alphaHW(E(i,1)) = accu(freqE2%(gvE2-gvEMu2)%(x-genoMu2))/
        accu(freqE2%(x-genoMu2)%(x-genoMu2));

      //Check for divide by zero
      if(!std::isfinite(alphaLe1)) alphaLe1=0;
      if(!std::isfinite(alphaHW(E(i,0)))) alphaHW(E(i,0))=0;
      if(!std::isfinite(alphaLe2)) alphaLe2=0;
      if(!std::isfinite(alphaHW(E(i,1)))) alphaHW(E(i,1))=0;

      //Breeding values
      arma::vec bv1, bv2, bvE1, bvE2;
      bv1 = (x-genoMu1)*alphaLe1; //Breeding values
      bvE1 = (x-genoMu1)*alphaHW(E(i,0)); //Random mating breeding value
      bv2 = (x-genoMu2)*alphaLe2; //Breeding values
      bvE2 = (x-genoMu2)*alphaHW(E(i,1)); //Random mating breeding value
      genicA(tid) += accu(freq1%bv1%bv1);
      genicA2(tid) += accu(freqE1%bvE1%bvE1);
      genicA(tid) += accu(freq2%bv2%bv2);
      genicA2(tid) += accu(freqE2%bvE2%bvE2);
      //Dominance deviation
      arma::vec dd1, dd2, ddE1, ddE2;
      if(hasD){
        dd1 = gv1-bv1-gvMu1; //Dominance deviations (lack of fit)
        ddE1 = gvE1-bvE1-gvEMu1; //Random mating dominance deviation
        dd2 = gv2-bv2-gvMu2; //Dominance deviations (lack of fit)
        ddE2 = gvE2-bvE2-gvEMu2; //Random mating dominance deviation
        genicD(tid) += accu(freq1%dd1%dd1);
        genicD2(tid) += accu(freqE1%ddE1%ddE1);
        genicD(tid) += accu(freq2%dd2%dd2);
        genicD2(tid) += accu(freqE2%ddE2%ddE2);
      }
    
      //Joint values (both loci)
      //Genetic value matrix
      arma::mat GV(ploidy+1,ploidy+1);
      //Breeding value matrix
      arma::mat BV(ploidy+1,ploidy+1);
      arma::mat BVE(ploidy+1,ploidy+1);
      //Dominance deviation matrix
      arma::mat DD(ploidy+1,ploidy+1);
      arma::mat DDE(ploidy+1,ploidy+1);
      //Epistasis matrix (lack of fit)
      arma::mat AA(ploidy+1,ploidy+1);
      arma::mat AANoLD(ploidy+1,ploidy+1);
      arma::mat AAE(ploidy+1,ploidy+1);
      for(arma::uword j=0; j<(ploidy+1); ++j){
        for(arma::uword k=0; k<(ploidy+1); ++k){
          BV(j,k) = bv1(j)+bv2(k);
          BVE(j,k) = bvE1(j)+bvE2(k);
          if(hasD){
            GV(j,k) = xa(j)*a(E(i,0)) + xa(k)*a(E(i,1)) +
              xd(j)*d(E(i,0)) + xd(k)*d(E(i,1)) + 
              xa(j)*xa(k)*E(i,2);
            DD(j,k) = dd1(j)+dd2(k);
            DDE(j,k) = ddE1(j)+ddE2(k);
          }else{
            GV(j,k) = xa(j)*a(E(i,0)) + xa(k)*a(E(i,1)) +
              xa(j)*xa(k)*E(i,2);
          }
        }
      }
      gvMu = accu(freq%GV);
      gvNoLDMu = accu(freqNoLD%GV);
      gvEMu = accu(freqNoLDE%GV);
      mu(tid) += gvMu;
      eMu(tid) += gvEMu;
      leMu(tid) += gvNoLDMu;
      if(hasD){
        AA = GV-BV-DD-gvMu;
        AANoLD = GV-BV-DD-gvNoLDMu;
        AAE = GV-BVE-DDE-gvEMu;
      }else{
        AA = GV-BV-gvMu;
        AANoLD = GV-BV-gvNoLDMu;
        AAE = GV-BVE-gvEMu;
      }
      genicAA(tid) += accu(freqNoLD%AANoLD%AANoLD);
      genicAA2(tid) += accu(freqNoLDE%AAE%AAE);

      if(locusTables){
        arma::uword t1 = loc1(i), t2 = loc2(i);
        tabFreq.col(t1) = freq1;
        tabFreq.col(t2) = freq2;
        tabFreqE.col(t1) = freqE1;
        tabFreqE.col(t2) = freqE2;
        tabBv.col(t1) = bv1;
        tabBv.col(t2) = bv2;
        tabBvE.col(t1) = bvE1;
        tabBvE.col(t2) = bvE2;
        if(hasD){
          tabDd.col(t1) = dd1;
          tabDd.col(t2) = dd2;
          tabDdE.col(t1) = ddE1;
          tabDdE.col(t2) = ddE2;
        }
        // Column major, so the dosage at the pair's first locus varies
        // fastest
        tabAA.col(i) = arma::vectorise(AANoLD);
        tabAAE.col(i) = arma::vectorise(AAE);
      }
    
      // Observed decomposition. Each locus's own effects are regressed on
      // its own dosage under its marginal frequencies, as for a trait
      // without epistasis.
      arma::uword l1 = loc1(i), l2 = loc2(i);
      arma::vec own1 = hasD ? arma::vec(aEff1+dEff1) : aEff1;
      arma::vec own2 = hasD ? arma::vec(aEff2+dEff2) : aEff2;
      double ownMu1 = accu(freq1%own1), ownMu2 = accu(freq2%own2);
      double alphaOwn1 = accu(freq1%(own1-ownMu1)%(x-genoMu1))/
        accu(freq1%(x-genoMu1)%(x-genoMu1));
      double alphaOwn2 = accu(freq2%(own2-ownMu2)%(x-genoMu2))/
        accu(freq2%(x-genoMu2)%(x-genoMu2));
      if(!std::isfinite(alphaOwn1)) alphaOwn1=0;
      if(!std::isfinite(alphaOwn2)) alphaOwn2=0;
      // The epistatic value of the pair is fitted to both dosages jointly
      // under the joint frequencies: on the dosages alone for its part of
      // the breeding values, and on the dosages and dominance codings for
      // its part of the dominance deviations. Column major order, so the
      // first locus's dosage varies fastest, as in freq.
      double xdMu1 = accu(freq1%xd), xdMu2 = accu(freq2%xd);
      arma::vec w = arma::vectorise(freq);
      arma::mat Z(nCell, 4);
      arma::vec yAA(nCell);
      for(arma::uword k=0; k<nGeno; ++k){
        for(arma::uword j=0; j<nGeno; ++j){
          arma::uword c = j+k*nGeno;
          Z(c,0) = x(j)-genoMu1;
          Z(c,1) = x(k)-genoMu2;
          Z(c,2) = xd(j)-xdMu1;
          Z(c,3) = xd(k)-xdMu2;
          yAA(c) = E(i,2)*xa(j)*xa(k);
        }
      }
      yAA -= accu(w%yAA);
      arma::vec gamma = wlsFit(Z, w, yAA, 2);
      arma::vec beta = wlsFit(Z, w, yAA, 4);
      alpha(l1) = alphaOwn1+gamma(0);
      alpha(l2) = alphaOwn2+gamma(1);
      // A locus's dominance deviation is its own effects less their
      // regression on its dosage, plus the part of the epistatic value the
      // dominance codings add to the joint fit, all written as effects on
      // the additive and dominance codings plus a constant, using
      // x = xa*P/2+P/2. Only the sum over the pair of the epistatic part is
      // uncorrelated with both dosages.
      double d1 = hasD ? d(l1) : 0.0, d2 = hasD ? d(l2) : 0.0;
      double extra1 = beta(0)-gamma(0), extra2 = beta(1)-gamma(1);
      locDdAdd(l1) = a(l1)+(extra1-alphaOwn1)*(dP/2.0);
      locDdAdd(l2) = a(l2)+(extra2-alphaOwn2)*(dP/2.0);
      locDdDom(l1) = d1+beta(2);
      locDdDom(l2) = d2+beta(3);
      locDdConst(l1) = -ownMu1+(extra1-alphaOwn1)*(dP/2.0-genoMu1)-
        beta(2)*xdMu1;
      locDdConst(l2) = -ownMu2+(extra2-alphaOwn2)*(dP/2.0-genoMu2)-
        beta(3)*xdMu2;
      locGenoMu(l1) = genoMu1;
      locGenoMu(l2) = genoMu2;
      pairGvMu(i) = gvMu;
    }
  }

  // Values for individuals. A locus's breeding value is alpha*(x-genoMu),
  // an effect of alpha*P/2 on the additive coding plus the constant
  // alpha*(P/2-genoMu), and its dominance deviation is a set of additive
  // and dominance effects plus a constant, worked out above.
  std::vector<GvEffects> effs;
  GvEffects eGvA, eBv, eGvD, eDd;
  eGvA.add = a;
  eGvA.hasA = true;
  eGvA.hasD = false;
  eGvA.hasGxe = false;
  eGvA.intercept = 0.0;
  eGvA.gxeInt = 0.0;
  eBv = eGvA;
  eBv.add = alpha*(dP/2.0);
  double bvConst = 0.0, ddConst = 0.0, pairMuSum = 0.0;
  for(arma::uword l=0; l<nLoci; ++l){
    bvConst += alpha(l)*(dP/2.0-locGenoMu(l));
    ddConst += locDdConst(l);
  }
  for(arma::uword i=0; i<nPair; ++i){
    pairMuSum += pairGvMu(i);
  }
  eBv.intercept = bvConst;
  effs.push_back(eGvA);
  effs.push_back(eBv);
  // gv_d is a genotypic value, so it exists only with dominance effects,
  // but under linkage disequilibrium the epistatic value can give dominance
  // deviations to any trait, so dd is always worked out and comes last
  if(hasD){
    eGvD = eGvA;
    eGvD.add.reset();
    eGvD.hasA = false;
    eGvD.dom = d;
    eGvD.hasD = true;
    effs.push_back(eGvD);
  }
  eDd = eGvA;
  eDd.add = locDdAdd;
  eDd.dom = locDdDom;
  eDd.hasD = true;
  eDd.intercept = ddConst;
  effs.push_back(eDd);
  std::vector<arma::field<arma::vec> > ind = gvGroupEffects(effs, g, pos,
                                                            nThreads);
  const arma::vec& gvA = ind[0](0);
  const arma::vec& bvInd = ind[1](0);
  const arma::vec& ddInd = ind[effs.size()-1](0);

  // Left empty unless requested
  Rcpp::List tables;
  if(locusTables){
    arma::uvec pair1 = loc1+1, pair2 = loc2+1; // C++ to R
    if(hasD){
      tables = Rcpp::List::create(Rcpp::Named("freq")=tabFreq,
                                  Rcpp::Named("freqE")=tabFreqE,
                                  Rcpp::Named("bv")=tabBv,
                                  Rcpp::Named("bvE")=tabBvE,
                                  Rcpp::Named("dd")=tabDd,
                                  Rcpp::Named("ddE")=tabDdE,
                                  Rcpp::Named("aa")=tabAA,
                                  Rcpp::Named("aaE")=tabAAE,
                                  Rcpp::Named("pair1")=pair1,
                                  Rcpp::Named("pair2")=pair2);
    }else{
      tables = Rcpp::List::create(Rcpp::Named("freq")=tabFreq,
                                  Rcpp::Named("freqE")=tabFreqE,
                                  Rcpp::Named("bv")=tabBv,
                                  Rcpp::Named("bvE")=tabBvE,
                                  Rcpp::Named("aa")=tabAA,
                                  Rcpp::Named("aaE")=tabAAE,
                                  Rcpp::Named("pair1")=pair1,
                                  Rcpp::Named("pair2")=pair2);
    }
  }

  if(hasD){
    const arma::vec& gvD = ind[2](0);
    return Rcpp::List::create(Rcpp::Named("gv")=gvA+gvD+gvAA+intercept,
                              Rcpp::Named("bv")=bvInd,
                              Rcpp::Named("dd")=ddInd,
                              Rcpp::Named("aa")=gvA+gvD+gvAA-bvInd-ddInd-pairMuSum,
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarD")=accu(genicD),
                              Rcpp::Named("genicVarAA")=accu(genicAA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarD2")=accu(genicD2),
                              Rcpp::Named("genicVarAA2")=accu(genicAA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("mu_LE")=accu(leMu)+intercept,
                              Rcpp::Named("gv_a")=gvA,
                              Rcpp::Named("gv_d")=gvD,
                              Rcpp::Named("gv_aa")=gvAA,
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW,
                              Rcpp::Named("tables")=tables);
  }else{
    return Rcpp::List::create(Rcpp::Named("gv")=gvA+gvAA+intercept,
                              Rcpp::Named("bv")=bvInd,
                              Rcpp::Named("dd")=ddInd,
                              Rcpp::Named("aa")=gvA+gvAA-bvInd-ddInd-pairMuSum,
                              Rcpp::Named("genicVarA")=accu(genicA),
                              Rcpp::Named("genicVarAA")=accu(genicAA),
                              Rcpp::Named("genicVarA2")=accu(genicA2),
                              Rcpp::Named("genicVarAA2")=accu(genicAA2),
                              Rcpp::Named("mu")=accu(mu)+intercept,
                              Rcpp::Named("mu_HWE")=accu(eMu)+intercept,
                              Rcpp::Named("mu_LE")=accu(leMu)+intercept,
                              Rcpp::Named("gv_a")=gvA,
                              Rcpp::Named("gv_aa")=gvAA,
                              Rcpp::Named("gv_mu")=intercept,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("alpha_HW")=alphaHW,
                              Rcpp::Named("tables")=tables);
  }
}

// Counts the genotypes at each locus: counts(x,i) is the number of
// individuals with dosage x at locus i. Loci sharing a byte are counted
// together, so each byte of each haplotype is read once per individual.
// Counts are integers, so they are exact whichever thread adds them.
arma::Mat<arma::uword> countGenotypes(const PackedGeno& g, const LociPos& pos,
                                      int nThreads){
  arma::uword nInd = g.nInd;
  arma::uword ploidy = g.ploidy;
  arma::uword nGeno = ploidy+1;
  arma::uword nLoci = pos.base.size();
  arma::Mat<arma::uword> counts(nGeno, nLoci, arma::fill::zeros);
  if(nLoci==0 || nInd==0){
    return counts;
  }

  // Loci in order of chromosome and byte, with the first locus of each
  // byte marked
  std::vector<arma::uword> cumBins(g.chrPtr.size()+1, 0);
  for(arma::uword c=0; c<g.chrPtr.size(); ++c){
    cumBins[c+1] = cumBins[c]+g.nBins[c];
  }
  std::vector<std::pair<arma::uword,arma::uword> > order(nLoci);
  for(arma::uword i=0; i<nLoci; ++i){
    order[i] = std::make_pair(cumBins[pos.chr(i)]+pos.byte(i), i);
  }
  std::sort(order.begin(), order.end());
  std::vector<arma::uword> keyStart;
  for(arma::uword i=0; i<nLoci; ++i){
    if(i==0 || order[i].first!=order[i-1].first){
      keyStart.push_back(i);
    }
  }
  arma::uword nKey = keyStart.size();
  keyStart.push_back(nLoci);

  // Each block owns the counters of its own loci
  arma::uword nBlocks = countBlocks(nKey);
  if(nBlocks < static_cast<arma::uword>(nThreads)){
    nThreads = static_cast<int>(nBlocks);
  }
  arma::uword* cnt = counts.memptr();
  if(ploidy==2){
    // A diploid's genotype at a locus follows from whether its two bits
    // differ and whether both are set. Those are worked out for a whole
    // byte at once, and spread[b] puts bit k of b into byte k of a 64 bit
    // word, so one addition counts all eight loci of a byte. The byte
    // counters are flushed every 255 individuals, before they can overflow.
    std::vector<uint64_t> spread(256);
    for(arma::uword b=0; b<256; ++b){
      uint64_t w = 0;
      for(arma::uword k=0; k<8; ++k){
        w |= uint64_t((b>>k) & 1u) << (8*k);
      }
      spread[b] = w;
    }
    const arma::uword kTile = 255;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for(arma::uword block=0; block<nBlocks; ++block){
      arma::uword kStart = blockStart(nKey, nBlocks, block);
      arma::uword kEnd = blockStart(nKey, nBlocks, block+1);
      for(arma::uword t0=0; t0<nInd; t0+=kTile){
        arma::uword t1 = std::min(t0+kTile, nInd);
        for(arma::uword k=kStart; k<kEnd; ++k){
          arma::uword first = order[keyStart[k]].second;
          const unsigned char* base = pos.base[first];
          arma::uword hs = pos.hapStride(first);
          arma::uword is = pos.indStride(first);
          uint64_t het = 0, hom = 0;
          for(arma::uword j=t0; j<t1; ++j){
            const unsigned char* p = base+is*j;
            unsigned char b0 = p[0], b1 = p[hs];
            het += spread[b0 ^ b1];
            hom += spread[b0 & b1];
          }
          for(arma::uword m=keyStart[k]; m<keyStart[k+1]; ++m){
            arma::uword i = order[m].second;
            arma::uword shift = 8*pos.bit(i);
            cnt[i*nGeno+1] += arma::uword((het>>shift) & 0xFFu);
            cnt[i*nGeno+2] += arma::uword((hom>>shift) & 0xFFu);
          }
        }
      }
    }
    for(arma::uword i=0; i<nLoci; ++i){
      cnt[i*nGeno] = nInd-cnt[i*nGeno+1]-cnt[i*nGeno+2];
    }
    return counts;
  }

  // Any other ploidy counts each locus's dosage directly
  const arma::uword kTile = 256;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword block=0; block<nBlocks; ++block){
    arma::uword kStart = blockStart(nKey, nBlocks, block);
    arma::uword kEnd = blockStart(nKey, nBlocks, block+1);
    std::vector<unsigned char> hapByte(ploidy);
    for(arma::uword t0=0; t0<nInd; t0+=kTile){
      arma::uword t1 = std::min(t0+kTile, nInd);
      for(arma::uword k=kStart; k<kEnd; ++k){
        arma::uword first = order[keyStart[k]].second;
        const unsigned char* base = pos.base[first];
        arma::uword hs = pos.hapStride(first);
        arma::uword is = pos.indStride(first);
        for(arma::uword j=t0; j<t1; ++j){
          const unsigned char* p = base+is*j;
          for(arma::uword h=0; h<ploidy; ++h){
            hapByte[h] = p[hs*h];
          }
          for(arma::uword m=keyStart[k]; m<keyStart[k+1]; ++m){
            arma::uword i = order[m].second;
            arma::uword bit = pos.bit(i);
            unsigned int dose = 0;
            for(arma::uword h=0; h<ploidy; ++h){
              dose += (unsigned int) ((hapByte[h]>>bit) & 1u);
            }
            ++cnt[i*nGeno+dose];
          }
        }
      }
    }
  }
  return counts;
}

// Calculates breeding values, dominance deviations and genic variances for
// traits without epistasis that all have exactly the loci in pos, as
// correlated traits added together do. Additive and dominance genetic
// variances are calculated from breeding values and dominance deviations.
//
// Genotypes are not decoded. The genotypes at each locus are counted from
// the packed bits once for all the traits, the statistics for each locus
// follow from those counts, and the values for individuals are then worked
// out the way genetic values are. Once a locus's genotype frequencies are
// known, an individual's breeding value at that locus is a function of its
// dosage with the same form as a trait's additive effects plus a constant,
// so gv_a, gv_d and bv for every trait are sets of effects on the shared
// loci and are worked out together by gvGroupEffects, which reads the
// genotypes once for all of them. Each trait's statistics for loci are the
// same whichever traits it is worked out with.
//
// With locusTables, the genotype frequencies and statistical values at each
// locus are also returned, so that genParam can form genic covariances
// between traits.
std::vector<Rcpp::List> calcGenParamGroup(const std::vector<Rcpp::S4>& traits,
                                          const PackedGeno& g,
                                          const LociPos& pos,
                                          int nThreads,
                                          bool locusTables=false){
  arma::uword nTrait = traits.size();
  arma::uword ploidy = g.ploidy;
  double dP = double(ploidy);
  arma::uword nLoci = pos.base.size();
  arma::vec x(ploidy+1); // Genotype dossage
  for(arma::uword i=0; i<x.n_elem; ++i)
    x(i) = double(i);
  arma::vec xa = (x-dP/2.0)*(2.0/dP);
  arma::vec xd = x%(dP-x)*(2.0/dP)*(2.0/dP);

  arma::Mat<arma::uword> counts = countGenotypes(g, pos, nThreads);

  // Sums over loci are accumulated by block of loci, so that their order
  // does not depend on the thread count
  arma::uword nBlocks = countBlocks(nLoci);
  int nThreadsLoci = nThreads;
  if(nBlocks < static_cast<arma::uword>(nThreadsLoci)){
    nThreadsLoci = static_cast<int>(nBlocks);
  }

  // Per trait results that do not depend on individuals
  std::vector<bool> hasD(nTrait);
  std::vector<double> intercept(nTrait), gvMuSum(nTrait);
  std::vector<arma::vec> genicA(nTrait), genicA2(nTrait), genicD(nTrait),
    genicD2(nTrait), mu(nTrait), eMu(nTrait), alpha(nTrait), alphaHW(nTrait);
  std::vector<GvEffects> effs;
  std::vector<arma::uword> gvAIdx(nTrait), bvIdx(nTrait), gvDIdx(nTrait);
  // Tables for genic covariances between traits, one column per locus.
  // Genotype frequencies belong to the loci, not the traits, so they are
  // written while the first trait is worked out and shared by the rest.
  // Each column is written by the one work item that owns the locus, and
  // none of these values feeds the sums.
  arma::uword nGeno = ploidy+1;
  arma::mat tabFreq, tabFreqE;
  std::vector<arma::mat> tabBv(nTrait), tabBvE(nTrait), tabDd(nTrait),
    tabDdE(nTrait);
  if(locusTables){
    tabFreq.set_size(nGeno, nLoci);
    tabFreqE.set_size(nGeno, nLoci);
  }

  for(arma::uword t=0; t<nTrait; ++t){
    const Rcpp::S4& trait = traits[t];
    hasD[t] = trait.hasSlot("domEff");
    arma::vec a = trait.slot("addEff");
    arma::vec d;
    if(hasD[t]){
      d = Rcpp::as<arma::vec>(trait.slot("domEff"));
    }
    if(a.n_elem!=nLoci){
      Rcpp::stop("number of loci does not match the number of effects");
    }
    intercept[t] = trait.slot("intercept");
    genicA[t].zeros(nBlocks); // No LD
    genicA2[t].zeros(nBlocks); // No LD and HWE
    genicD[t].zeros(nBlocks); // No LD
    genicD2[t].zeros(nBlocks); // No LD and HWE
    mu[t].zeros(nBlocks); // Observed mean
    eMu[t].zeros(nBlocks); // Expected mean with HWE
    alpha[t].set_size(nLoci);
    alphaHW[t].set_size(nLoci);
    arma::vec locusGenoMu(nLoci), locusGvMu(nLoci);
    bool hasDt = hasD[t];
    arma::vec& genicAt = genicA[t];
    arma::vec& genicA2t = genicA2[t];
    arma::vec& genicDt = genicD[t];
    arma::vec& genicD2t = genicD2[t];
    arma::vec& mut = mu[t];
    arma::vec& eMut = eMu[t];
    arma::vec& alphat = alpha[t];
    arma::vec& alphaHWt = alphaHW[t];
    if(locusTables){
      tabBv[t].set_size(nGeno, nLoci);
      tabBvE[t].set_size(nGeno, nLoci);
      if(hasDt){
        tabDd[t].set_size(nGeno, nLoci);
        tabDdE[t].set_size(nGeno, nLoci);
      }
    }
    arma::mat& tabBvt = tabBv[t];
    arma::mat& tabBvEt = tabBvE[t];
    arma::mat& tabDdt = tabDd[t];
    arma::mat& tabDdEt = tabDdE[t];
    bool writeFreq = locusTables && t==0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreadsLoci)
#endif
    for(arma::uword tid=0; tid<nBlocks; ++tid){
      arma::uword itemStart = blockStart(nLoci, nBlocks, tid);
      arma::uword itemEnd = blockStart(nLoci, nBlocks, tid+1);
      arma::vec freq(ploidy+1), freqE(ploidy+1); // Genotype frequencies, observed and HWE
      arma::vec aEff(ploidy+1), dEff(ploidy+1), gv(ploidy+1); // Genetic values, additive and dominance
      arma::vec bv(ploidy+1), dd(ploidy+1); // Statistical values, additive and dominance
      arma::vec bvE(ploidy+1), ddE(ploidy+1); //Expected for random mating
      for(arma::uword i=itemStart; i<itemEnd; ++i){
        double gvMu, gvEMu, genoMu, p, q, dK;

        // Genotype frequencies
        for(arma::uword k=0; k<(ploidy+1); ++k){
          freq(k) = double(counts(k,i));
        }
        freq = freq/accu(freq);
        genoMu = accu(freq%x);
        p = genoMu/dP;
        q = 1-p;

        // Expected genotype frequencies
        freqE.zeros();
        for(arma::uword k=0; k<(ploidy+1); ++k){
          dK = double(k);
          freqE(k) = choose(dP,dK)*std::pow(p,dK)*std::pow(q,dP-dK);
        }

        // Set genetic values. The two parts are worked out separately, as
        // before, so that the compiler cannot fuse them into one rounding.
        aEff = xa*a(i);
        if(hasDt){
          dEff = xd*d(i);
          gv = aEff+dEff;
        }else{
          gv = aEff;
        }

        // Mean genetic values
        gvMu = accu(freq%gv);
        gvEMu =  accu(freqE%gv);
        mut(tid) += gvMu;
        eMut(tid) += gvEMu;

        // Average effect
        alphat(i) = accu(freq%(gv-gvMu)%(x-genoMu))/
          accu(freq%(x-genoMu)%(x-genoMu));
        alphaHWt(i) = accu(freqE%(gv-gvEMu)%(x-genoMu))/
          accu(freqE%(x-genoMu)%(x-genoMu));

        // Check for divide by zero
        if(!std::isfinite(alphat(i))) alphat(i)=0;
        if(!std::isfinite(alphaHWt(i))) alphaHWt(i)=0;

        // Set additive genic variances
        bv = (x-genoMu)*alphat(i); //Breeding values
        bvE = (x-genoMu)*alphaHWt(i); //Random mating breeding value
        genicAt(tid) += accu(freq%bv%bv);
        genicA2t(tid) += accu(freqE%bvE%bvE);

        // Set dominance genic variances
        if(hasDt){
          dd = gv-bv-gvMu; //Dominance deviations (lack of fit)
          ddE = gv-bvE-gvEMu; //Random mating dominance deviation
          genicDt(tid) += accu(freq%dd%dd);
          genicD2t(tid) += accu(freqE%ddE%ddE);
        }

        locusGenoMu(i) = genoMu;
        locusGvMu(i) = gvMu;

        if(locusTables){
          tabBvt.col(i) = bv;
          tabBvEt.col(i) = bvE;
          if(hasDt){
            tabDdt.col(i) = dd;
            tabDdEt.col(i) = ddE;
          }
          if(writeFreq){
            tabFreq.col(i) = freq;
            tabFreqE.col(i) = freqE;
          }
        }
      }
    }

    // Values for individuals. The breeding value at a locus is
    // alpha*(x-genoMu), which on the additive coding is an effect of
    // alpha*P/2 plus the constant alpha*(P/2-genoMu). The constants are
    // summed over loci into the intercept. The dominance deviation at a
    // locus is gv-bv-gvMu, so an individual's dominance deviation is
    // gv_a+gv_d-bv less the sum of gvMu, and does not need a pass of its own.
    GvEffects eGvA, eBv, eGvD;
    eGvA.add = a;
    eGvA.hasA = true;
    eGvA.hasD = false;
    eGvA.hasGxe = false;
    eGvA.intercept = 0.0;
    eGvA.gxeInt = 0.0;
    eBv = eGvA;
    eBv.add = alphat*(dP/2.0);
    double bvConst = 0.0;
    gvMuSum[t] = 0.0;
    for(arma::uword i=0; i<nLoci; ++i){
      bvConst += alphat(i)*(dP/2.0-locusGenoMu(i));
      gvMuSum[t] += locusGvMu(i);
    }
    eBv.intercept = bvConst;
    gvAIdx[t] = effs.size();
    effs.push_back(eGvA);
    bvIdx[t] = effs.size();
    effs.push_back(eBv);
    if(hasDt){
      eGvD = eGvA;
      eGvD.add.reset();
      eGvD.hasA = false;
      eGvD.dom = d;
      eGvD.hasD = true;
      gvDIdx[t] = effs.size();
      effs.push_back(eGvD);
    }
  }

  std::vector<arma::field<arma::vec> > ind = gvGroupEffects(effs, g, pos,
                                                            nThreads);

  std::vector<Rcpp::List> output(nTrait);
  for(arma::uword t=0; t<nTrait; ++t){
    const arma::vec& gvA = ind[gvAIdx[t]](0);
    const arma::vec& bvInd = ind[bvIdx[t]](0);
    // Left empty unless requested
    Rcpp::List tables;
    if(locusTables){
      if(hasD[t]){
        tables = Rcpp::List::create(Rcpp::Named("freq")=tabFreq,
                                    Rcpp::Named("freqE")=tabFreqE,
                                    Rcpp::Named("bv")=tabBv[t],
                                    Rcpp::Named("bvE")=tabBvE[t],
                                    Rcpp::Named("dd")=tabDd[t],
                                    Rcpp::Named("ddE")=tabDdE[t]);
      }else{
        tables = Rcpp::List::create(Rcpp::Named("freq")=tabFreq,
                                    Rcpp::Named("freqE")=tabFreqE,
                                    Rcpp::Named("bv")=tabBv[t],
                                    Rcpp::Named("bvE")=tabBvE[t]);
      }
    }
    // Without epistasis the genetic value is a sum of functions of single
    // loci, so its mean does not depend on linkage disequilibrium and the
    // mean with no LD, mu_LE, is the observed mean
    if(hasD[t]){
      const arma::vec& gvD = ind[gvDIdx[t]](0);
      output[t] = Rcpp::List::create(Rcpp::Named("gv")=gvA+gvD+intercept[t],
                                     Rcpp::Named("bv")=bvInd,
                                     Rcpp::Named("dd")=gvA+gvD-bvInd-gvMuSum[t],
                                     Rcpp::Named("genicVarA")=accu(genicA[t]),
                                     Rcpp::Named("genicVarD")=accu(genicD[t]),
                                     Rcpp::Named("genicVarA2")=accu(genicA2[t]),
                                     Rcpp::Named("genicVarD2")=accu(genicD2[t]),
                                     Rcpp::Named("mu")=accu(mu[t])+intercept[t],
                                     Rcpp::Named("mu_HWE")=accu(eMu[t])+intercept[t],
                                     Rcpp::Named("mu_LE")=accu(mu[t])+intercept[t],
                                     Rcpp::Named("gv_a")=gvA,
                                     Rcpp::Named("gv_d")=gvD,
                                     Rcpp::Named("gv_mu")=intercept[t],
                                     Rcpp::Named("alpha")=alpha[t],
                                     Rcpp::Named("alpha_HW")=alphaHW[t],
                                     Rcpp::Named("tables")=tables);
    }else{
      output[t] = Rcpp::List::create(Rcpp::Named("gv")=gvA+intercept[t],
                                     Rcpp::Named("bv")=bvInd,
                                     Rcpp::Named("genicVarA")=accu(genicA[t]),
                                     Rcpp::Named("genicVarA2")=accu(genicA2[t]),
                                     Rcpp::Named("mu")=accu(mu[t])+intercept[t],
                                     Rcpp::Named("mu_HWE")=accu(eMu[t])+intercept[t],
                                     Rcpp::Named("mu_LE")=accu(mu[t])+intercept[t],
                                     Rcpp::Named("gv_a")=gvA,
                                     Rcpp::Named("gv_mu")=intercept[t],
                                     Rcpp::Named("alpha")=alpha[t],
                                     Rcpp::Named("alpha_HW")=alphaHW[t],
                                     Rcpp::Named("tables")=tables);
    }
  }
  return output;
}

// Calculates breeding values, dominance deviations and genic
// variances for one trait. Additive and dominance genetic variances are
// calculated from breeding values and dominance deviations.
// [[Rcpp::export]]
Rcpp::List calcGenParam(const Rcpp::S4& trait,
                        const Rcpp::S4& pop,
                        int nThreads){
  if(trait.hasSlot("epiEff")){
    return calcGenParamE(trait, pop, nThreads);
  }
  PackedGeno g(pop);
  LociPos pos = findLoci(g, Rcpp::as<arma::Col<int> >(trait.slot("lociPerChr")),
                         Rcpp::as<arma::uvec>(trait.slot("lociLoc")));
  return calcGenParamGroup(std::vector<Rcpp::S4>(1, trait), g, pos,
                           nThreads)[0];
}

// calcGenParam for every trait in traitList, returned in a list in the same
// order. Traits without epistasis that have exactly the same loci, as
// correlated traits added together do, are worked out together, so their
// genotypes are counted and read once for all of them. With locusTables,
// each trait's result also holds the tables genParam uses to form genic
// covariances between traits.
// [[Rcpp::export]]
Rcpp::List calcGenParamMulti(const Rcpp::List& traitList,
                             const Rcpp::S4& pop,
                             int nThreads,
                             bool locusTables=false){
  R_xlen_t nTraits = traitList.size();
  Rcpp::List output(nTraits);
  PackedGeno g(pop);
  std::vector<std::vector<R_xlen_t> > groups;
  std::vector<arma::Col<int> > groupPerChr;
  std::vector<arma::uvec> groupLoc;
  for(R_xlen_t i=0; i<nTraits; ++i){
    Rcpp::S4 trait = traitList[i];
    if(trait.hasSlot("epiEff")){
      output[i] = calcGenParamE(trait, pop, nThreads, locusTables);
      continue;
    }
    arma::Col<int> perChr = Rcpp::as<arma::Col<int> >(trait.slot("lociPerChr"));
    arma::uvec loc = Rcpp::as<arma::uvec>(trait.slot("lociLoc"));
    bool placed = false;
    for(arma::uword k=0; k<groups.size() && !placed; ++k){
      if(groupPerChr[k].n_elem==perChr.n_elem &&
         arma::all(groupPerChr[k]==perChr) &&
         groupLoc[k].n_elem==loc.n_elem &&
         arma::all(groupLoc[k]==loc)){
        groups[k].push_back(i);
        placed = true;
      }
    }
    if(!placed){
      groups.push_back(std::vector<R_xlen_t>(1, i));
      groupPerChr.push_back(perChr);
      groupLoc.push_back(loc);
    }
  }
  for(arma::uword k=0; k<groups.size(); ++k){
    std::vector<Rcpp::S4> members;
    for(arma::uword t=0; t<groups[k].size(); ++t){
      members.push_back(traitList[groups[k][t]]);
    }
    LociPos pos = findLoci(g, groupPerChr[k], groupLoc[k]);
    std::vector<Rcpp::List> res = calcGenParamGroup(members, g, pos, nThreads,
                                                    locusTables);
    for(arma::uword t=0; t<groups[k].size(); ++t){
      output[groups[k][t]] = res[t];
    }
  }
  return output;
}
