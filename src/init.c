#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _AlphaSimR_calcChrFreq(SEXP);
extern SEXP _AlphaSimR_calcCoef(SEXP, SEXP);
extern SEXP _AlphaSimR_calcGenParam(SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_GCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_MV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_SCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_convToImat(SEXP);
extern SEXP _AlphaSimR_createDH2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_cross2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_gebvGCA(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_gebvRR(SEXP, SEXP);
extern SEXP _AlphaSimR_gegvGCA(SEXP, SEXP);
extern SEXP _AlphaSimR_gegvRRD(SEXP, SEXP);
extern SEXP _AlphaSimR_gegvSCA(SEXP, SEXP);
extern SEXP _AlphaSimR_getDomGeno(SEXP);
extern SEXP _AlphaSimR_getGeno(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getGv(SEXP, SEXP);
extern SEXP _AlphaSimR_getHaplo(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getHybridGv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getIbdHaplo(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getIbdRecHist(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getOneHaplo(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_MaCS(SEXP, SEXP);
extern SEXP _AlphaSimR_mergeGeno(SEXP, SEXP);
extern SEXP _AlphaSimR_mergeMultGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_mergeMultIntMat(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_packHaplo(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_popVar(SEXP);
extern SEXP _AlphaSimR_sampAllComb(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_sampHalfDialComb(SEXP, SEXP);
extern SEXP _AlphaSimR_tuneTraitA(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_tuneTraitAD(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeASGenotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeASHaplotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeOneHaplo(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writePlinkPed(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_AlphaSimR_calcChrFreq",       (DL_FUNC) &_AlphaSimR_calcChrFreq,        1},
    {"_AlphaSimR_calcCoef",          (DL_FUNC) &_AlphaSimR_calcCoef,           2},
    {"_AlphaSimR_calcGenParam",      (DL_FUNC) &_AlphaSimR_calcGenParam,       2},
    {"_AlphaSimR_callRRBLUP",        (DL_FUNC) &_AlphaSimR_callRRBLUP,         6},
    {"_AlphaSimR_callRRBLUP_D",      (DL_FUNC) &_AlphaSimR_callRRBLUP_D,       8},
    {"_AlphaSimR_callRRBLUP_GCA",    (DL_FUNC) &_AlphaSimR_callRRBLUP_GCA,     8},
    {"_AlphaSimR_callRRBLUP_MV",     (DL_FUNC) &_AlphaSimR_callRRBLUP_MV,      7},
    {"_AlphaSimR_callRRBLUP_SCA",    (DL_FUNC) &_AlphaSimR_callRRBLUP_SCA,     9},
    {"_AlphaSimR_callRRBLUP2",       (DL_FUNC) &_AlphaSimR_callRRBLUP2,       11},
    {"_AlphaSimR_convToImat",        (DL_FUNC) &_AlphaSimR_convToImat,         1},
    {"_AlphaSimR_createDH2",         (DL_FUNC) &_AlphaSimR_createDH2,          4},
    {"_AlphaSimR_cross2",            (DL_FUNC) &_AlphaSimR_cross2,             7},
    {"_AlphaSimR_gebvGCA",           (DL_FUNC) &_AlphaSimR_gebvGCA,            3},
    {"_AlphaSimR_gebvRR",            (DL_FUNC) &_AlphaSimR_gebvRR,             2},
    {"_AlphaSimR_gegvGCA",           (DL_FUNC) &_AlphaSimR_gegvGCA,            2},
    {"_AlphaSimR_gegvRRD",           (DL_FUNC) &_AlphaSimR_gegvRRD,            2},
    {"_AlphaSimR_gegvSCA",           (DL_FUNC) &_AlphaSimR_gegvSCA,            2},
    {"_AlphaSimR_getDomGeno",        (DL_FUNC) &_AlphaSimR_getDomGeno,         1},
    {"_AlphaSimR_getGeno",           (DL_FUNC) &_AlphaSimR_getGeno,            3},
    {"_AlphaSimR_getGv",             (DL_FUNC) &_AlphaSimR_getGv,              2},
    {"_AlphaSimR_getHaplo",          (DL_FUNC) &_AlphaSimR_getHaplo,           3},
    {"_AlphaSimR_getHybridGv",       (DL_FUNC) &_AlphaSimR_getHybridGv,        5},
    {"_AlphaSimR_getIbdHaplo",       (DL_FUNC) &_AlphaSimR_getIbdHaplo,        3},
    {"_AlphaSimR_getIbdRecHist",     (DL_FUNC) &_AlphaSimR_getIbdRecHist,      4},
    {"_AlphaSimR_getOneHaplo",       (DL_FUNC) &_AlphaSimR_getOneHaplo,        4},
    {"_AlphaSimR_MaCS",              (DL_FUNC) &_AlphaSimR_MaCS,               2},
    {"_AlphaSimR_mergeGeno",         (DL_FUNC) &_AlphaSimR_mergeGeno,          2},
    {"_AlphaSimR_mergeMultGeno",     (DL_FUNC) &_AlphaSimR_mergeMultGeno,      4},
    {"_AlphaSimR_mergeMultIntMat",   (DL_FUNC) &_AlphaSimR_mergeMultIntMat,    3},
    {"_AlphaSimR_packHaplo",         (DL_FUNC) &_AlphaSimR_packHaplo,          3},
    {"_AlphaSimR_popVar",            (DL_FUNC) &_AlphaSimR_popVar,             1},
    {"_AlphaSimR_sampAllComb",       (DL_FUNC) &_AlphaSimR_sampAllComb,        3},
    {"_AlphaSimR_sampHalfDialComb",  (DL_FUNC) &_AlphaSimR_sampHalfDialComb,   2},
    {"_AlphaSimR_tuneTraitA",        (DL_FUNC) &_AlphaSimR_tuneTraitA,         3},
    {"_AlphaSimR_tuneTraitAD",       (DL_FUNC) &_AlphaSimR_tuneTraitAD,        5},
    {"_AlphaSimR_writeASGenotypes",  (DL_FUNC) &_AlphaSimR_writeASGenotypes,   7},
    {"_AlphaSimR_writeASHaplotypes", (DL_FUNC) &_AlphaSimR_writeASHaplotypes,  7},
    {"_AlphaSimR_writeGeno",         (DL_FUNC) &_AlphaSimR_writeGeno,          4},
    {"_AlphaSimR_writeOneHaplo",     (DL_FUNC) &_AlphaSimR_writeOneHaplo,      5},
    {"_AlphaSimR_writePlinkPed",     (DL_FUNC) &_AlphaSimR_writePlinkPed,      5},
    {NULL, NULL, 0}
};

void R_init_AlphaSimR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
