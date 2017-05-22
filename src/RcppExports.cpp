// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// AlphaFormatter
int AlphaFormatter();
RcppExport SEXP AlphaSimR_AlphaFormatter() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(AlphaFormatter());
    return rcpp_result_gen;
END_RCPP
}
// getGeno
arma::Mat<unsigned char> getGeno(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc);
RcppExport SEXP AlphaSimR_getGeno(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    rcpp_result_gen = Rcpp::wrap(getGeno(geno, lociPerChr, lociLoc));
    return rcpp_result_gen;
END_RCPP
}
// getDomGeno
arma::imat getDomGeno(const arma::Mat<unsigned char>& geno);
RcppExport SEXP AlphaSimR_getDomGeno(SEXP genoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<unsigned char>& >::type geno(genoSEXP);
    rcpp_result_gen = Rcpp::wrap(getDomGeno(geno));
    return rcpp_result_gen;
END_RCPP
}
// getHaplo
arma::Mat<unsigned char> getHaplo(const arma::field<arma::Cube<unsigned char> >& geno, const arma::ivec& lociPerChr, arma::uvec lociLoc);
RcppExport SEXP AlphaSimR_getHaplo(SEXP genoSEXP, SEXP lociPerChrSEXP, SEXP lociLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type lociPerChr(lociPerChrSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type lociLoc(lociLocSEXP);
    rcpp_result_gen = Rcpp::wrap(getHaplo(geno, lociPerChr, lociLoc));
    return rcpp_result_gen;
END_RCPP
}
// getGvA
arma::vec getGvA(const Rcpp::S4& trait, const Rcpp::S4& pop);
RcppExport SEXP AlphaSimR_getGvA(SEXP traitSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(getGvA(trait, pop));
    return rcpp_result_gen;
END_RCPP
}
// getGvAG
arma::vec getGvAG(const Rcpp::S4& trait, const Rcpp::S4& pop, double z);
RcppExport SEXP AlphaSimR_getGvAG(SEXP traitSEXP, SEXP popSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(getGvAG(trait, pop, z));
    return rcpp_result_gen;
END_RCPP
}
// getGvAD
arma::vec getGvAD(const Rcpp::S4& trait, const Rcpp::S4& pop);
RcppExport SEXP AlphaSimR_getGvAD(SEXP traitSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(getGvAD(trait, pop));
    return rcpp_result_gen;
END_RCPP
}
// getGvADG
arma::vec getGvADG(const Rcpp::S4& trait, const Rcpp::S4& pop, double z);
RcppExport SEXP AlphaSimR_getGvADG(SEXP traitSEXP, SEXP popSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(getGvADG(trait, pop, z));
    return rcpp_result_gen;
END_RCPP
}
// calcGenParam
Rcpp::List calcGenParam(const Rcpp::S4& trait, const Rcpp::S4& pop);
RcppExport SEXP AlphaSimR_calcGenParam(SEXP traitSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(calcGenParam(trait, pop));
    return rcpp_result_gen;
END_RCPP
}
// getHybridGvA
arma::vec getHybridGvA(const Rcpp::S4& trait, const Rcpp::S4& fPop, arma::uvec& fPar, const Rcpp::S4& mPop, arma::uvec& mPar);
RcppExport SEXP AlphaSimR_getHybridGvA(SEXP traitSEXP, SEXP fPopSEXP, SEXP fParSEXP, SEXP mPopSEXP, SEXP mParSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type fPop(fPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type fPar(fParSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type mPop(mPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type mPar(mParSEXP);
    rcpp_result_gen = Rcpp::wrap(getHybridGvA(trait, fPop, fPar, mPop, mPar));
    return rcpp_result_gen;
END_RCPP
}
// getHybridGvAG
arma::vec getHybridGvAG(const Rcpp::S4& trait, const Rcpp::S4& fPop, arma::uvec& fPar, const Rcpp::S4& mPop, arma::uvec& mPar, double z);
RcppExport SEXP AlphaSimR_getHybridGvAG(SEXP traitSEXP, SEXP fPopSEXP, SEXP fParSEXP, SEXP mPopSEXP, SEXP mParSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type fPop(fPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type fPar(fParSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type mPop(mPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type mPar(mParSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(getHybridGvAG(trait, fPop, fPar, mPop, mPar, z));
    return rcpp_result_gen;
END_RCPP
}
// getHybridGvAD
arma::vec getHybridGvAD(const Rcpp::S4& trait, const Rcpp::S4& fPop, arma::uvec& fPar, const Rcpp::S4& mPop, arma::uvec& mPar);
RcppExport SEXP AlphaSimR_getHybridGvAD(SEXP traitSEXP, SEXP fPopSEXP, SEXP fParSEXP, SEXP mPopSEXP, SEXP mParSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type fPop(fPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type fPar(fParSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type mPop(mPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type mPar(mParSEXP);
    rcpp_result_gen = Rcpp::wrap(getHybridGvAD(trait, fPop, fPar, mPop, mPar));
    return rcpp_result_gen;
END_RCPP
}
// getHybridGvADG
arma::vec getHybridGvADG(const Rcpp::S4& trait, const Rcpp::S4& fPop, arma::uvec& fPar, const Rcpp::S4& mPop, arma::uvec& mPar, double z);
RcppExport SEXP AlphaSimR_getHybridGvADG(SEXP traitSEXP, SEXP fPopSEXP, SEXP fParSEXP, SEXP mPopSEXP, SEXP mParSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type trait(traitSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type fPop(fPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type fPar(fParSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type mPop(mPopSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type mPar(mParSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(getHybridGvADG(trait, fPop, fPar, mPop, mPar, z));
    return rcpp_result_gen;
END_RCPP
}
// cross2
arma::field<arma::Cube<unsigned char> > cross2(const arma::field<arma::Cube<unsigned char> >& fGeno, arma::uvec fPar, const arma::field<arma::Cube<unsigned char> >& mGeno, arma::uvec mPar, const arma::field<arma::vec>& genMaps);
RcppExport SEXP AlphaSimR_cross2(SEXP fGenoSEXP, SEXP fParSEXP, SEXP mGenoSEXP, SEXP mParSEXP, SEXP genMapsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type fGeno(fGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fPar(fParSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type mGeno(mGenoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mPar(mParSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::vec>& >::type genMaps(genMapsSEXP);
    rcpp_result_gen = Rcpp::wrap(cross2(fGeno, fPar, mGeno, mPar, genMaps));
    return rcpp_result_gen;
END_RCPP
}
// createDH2
arma::field<arma::Cube<unsigned char> > createDH2(const arma::field<arma::Cube<unsigned char> >& geno, int nDH, const arma::field<arma::vec>& genMaps);
RcppExport SEXP AlphaSimR_createDH2(SEXP genoSEXP, SEXP nDHSEXP, SEXP genMapsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< int >::type nDH(nDHSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::vec>& >::type genMaps(genMapsSEXP);
    rcpp_result_gen = Rcpp::wrap(createDH2(geno, nDH, genMaps));
    return rcpp_result_gen;
END_RCPP
}
// popVar
arma::mat popVar(const arma::mat& X);
RcppExport SEXP AlphaSimR_popVar(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(popVar(X));
    return rcpp_result_gen;
END_RCPP
}
// mergeGeno
arma::field<arma::Cube<unsigned char> > mergeGeno(const arma::field<arma::Cube<unsigned char> >& x, const arma::field<arma::Cube<unsigned char> >& y);
RcppExport SEXP AlphaSimR_mergeGeno(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::Cube<unsigned char> >& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mergeGeno(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calcChrMinorFreq
arma::vec calcChrMinorFreq(const arma::Cube<unsigned char>& geno, int ploidy);
RcppExport SEXP AlphaSimR_calcChrMinorFreq(SEXP genoSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Cube<unsigned char>& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(calcChrMinorFreq(geno, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// convToImat
arma::imat convToImat(const arma::Mat<unsigned char>& X);
RcppExport SEXP AlphaSimR_convToImat(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<unsigned char>& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(convToImat(X));
    return rcpp_result_gen;
END_RCPP
}
// solveUVM
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X, const arma::mat& Z, const arma::mat& K);
RcppExport SEXP AlphaSimR_solveUVM(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(solveUVM(y, X, Z, K));
    return rcpp_result_gen;
END_RCPP
}
// solveMVM
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X, const arma::mat& Z, const arma::mat& K, double tol);
RcppExport SEXP AlphaSimR_solveMVM(SEXP YSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP KSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMVM(Y, X, Z, K, tol));
    return rcpp_result_gen;
END_RCPP
}
// objWeights
Rcpp::NumericVector objWeights(Rcpp::NumericVector x, SEXP ptrData);
RcppExport SEXP AlphaSimR_objWeights(SEXP xSEXP, SEXP ptrDataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ptrData(ptrDataSEXP);
    rcpp_result_gen = Rcpp::wrap(objWeights(x, ptrData));
    return rcpp_result_gen;
END_RCPP
}
// solveMKM
Rcpp::List solveMKM(arma::mat& y, arma::mat& X, arma::field<arma::mat>& Zlist, arma::field<arma::mat>& Klist);
RcppExport SEXP AlphaSimR_solveMKM(SEXP ySEXP, SEXP XSEXP, SEXP ZlistSEXP, SEXP KlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type Zlist(ZlistSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat>& >::type Klist(KlistSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMKM(y, X, Zlist, Klist));
    return rcpp_result_gen;
END_RCPP
}
// callGK
Rcpp::List callGK(arma::mat y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker, double maxTheta, int maxIter, bool writeForPred);
RcppExport SEXP AlphaSimR_callGK(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP, SEXP maxThetaSEXP, SEXP maxIterSEXP, SEXP writeForPredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    Rcpp::traits::input_parameter< double >::type maxTheta(maxThetaSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type writeForPred(writeForPredSEXP);
    rcpp_result_gen = Rcpp::wrap(callGK(y, x, reps, genoTrain, nMarker, maxTheta, maxIter, writeForPred));
    return rcpp_result_gen;
END_RCPP
}
// callPredGK
arma::mat callPredGK(const Rcpp::DataFrame& genoPred);
RcppExport SEXP AlphaSimR_callPredGK(SEXP genoPredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type genoPred(genoPredSEXP);
    rcpp_result_gen = Rcpp::wrap(callPredGK(genoPred));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP
Rcpp::List callRRBLUP(arma::mat y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker);
RcppExport SEXP AlphaSimR_callRRBLUP(SEXP ySEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP(y, x, reps, genoTrain, nMarker));
    return rcpp_result_gen;
END_RCPP
}
// callRRBLUP_MV
Rcpp::List callRRBLUP_MV(arma::mat Y, arma::uvec x, arma::vec reps, std::string genoTrain, int nMarker);
RcppExport SEXP AlphaSimR_callRRBLUP_MV(SEXP YSEXP, SEXP xSEXP, SEXP repsSEXP, SEXP genoTrainSEXP, SEXP nMarkerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< std::string >::type genoTrain(genoTrainSEXP);
    Rcpp::traits::input_parameter< int >::type nMarker(nMarkerSEXP);
    rcpp_result_gen = Rcpp::wrap(callRRBLUP_MV(Y, x, reps, genoTrain, nMarker));
    return rcpp_result_gen;
END_RCPP
}
// readAF
arma::Cube<unsigned char> readAF(int nInd, int segSites, int ploidy, arma::uvec keep, bool inbred);
RcppExport SEXP AlphaSimR_readAF(SEXP nIndSEXP, SEXP segSitesSEXP, SEXP ploidySEXP, SEXP keepSEXP, SEXP inbredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nInd(nIndSEXP);
    Rcpp::traits::input_parameter< int >::type segSites(segSitesSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< bool >::type inbred(inbredSEXP);
    rcpp_result_gen = Rcpp::wrap(readAF(nInd, segSites, ploidy, keep, inbred));
    return rcpp_result_gen;
END_RCPP
}
// tuneTraitA
Rcpp::List tuneTraitA(arma::Mat<unsigned char>& geno, arma::vec& addEff, double varG);
RcppExport SEXP AlphaSimR_tuneTraitA(SEXP genoSEXP, SEXP addEffSEXP, SEXP varGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned char>& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type addEff(addEffSEXP);
    Rcpp::traits::input_parameter< double >::type varG(varGSEXP);
    rcpp_result_gen = Rcpp::wrap(tuneTraitA(geno, addEff, varG));
    return rcpp_result_gen;
END_RCPP
}
// tuneTraitAD
Rcpp::List tuneTraitAD(arma::Mat<unsigned char>& geno, arma::vec& addEff, arma::vec& domEff, double varG);
RcppExport SEXP AlphaSimR_tuneTraitAD(SEXP genoSEXP, SEXP addEffSEXP, SEXP domEffSEXP, SEXP varGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned char>& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type addEff(addEffSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type domEff(domEffSEXP);
    Rcpp::traits::input_parameter< double >::type varG(varGSEXP);
    rcpp_result_gen = Rcpp::wrap(tuneTraitAD(geno, addEff, domEff, varG));
    return rcpp_result_gen;
END_RCPP
}
