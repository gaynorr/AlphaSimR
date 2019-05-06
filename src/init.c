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
extern SEXP _AlphaSimR_calcGenParam(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callFastRRBLUP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_D2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_GCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_GCA2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_MV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP_SCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_callRRBLUP2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_createDH2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_cross(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getGv(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getHaplo(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getHybridGv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getIbdHaplo(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getIbdRecHist(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getMaternalGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getNumThreads();
extern SEXP _AlphaSimR_getOneHaplo(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_getPaternalGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_MaCS(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_mergeGeno(SEXP, SEXP);
extern SEXP _AlphaSimR_mergeMultGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_mergeMultIntMat(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_packHaplo(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_popVar(SEXP);
extern SEXP _AlphaSimR_sampAllComb(SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_sampHalfDialComb(SEXP, SEXP);
extern SEXP _AlphaSimR_writeASGenotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeASHaplotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeGeno(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writeOneHaplo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _AlphaSimR_writePlinkPed(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_AlphaSimR_calcChrFreq",       (DL_FUNC) &_AlphaSimR_calcChrFreq,        1},
    {"_AlphaSimR_calcCoef",          (DL_FUNC) &_AlphaSimR_calcCoef,           2},
    {"_AlphaSimR_calcGenParam",      (DL_FUNC) &_AlphaSimR_calcGenParam,       3},
    {"_AlphaSimR_callFastRRBLUP",    (DL_FUNC) &_AlphaSimR_callFastRRBLUP,     8},
    {"_AlphaSimR_callRRBLUP",        (DL_FUNC) &_AlphaSimR_callRRBLUP,         7},
    {"_AlphaSimR_callRRBLUP_D",      (DL_FUNC) &_AlphaSimR_callRRBLUP_D,       8},
    {"_AlphaSimR_callRRBLUP_D2",     (DL_FUNC) &_AlphaSimR_callRRBLUP_D2,     13},
    {"_AlphaSimR_callRRBLUP_GCA",    (DL_FUNC) &_AlphaSimR_callRRBLUP_GCA,     8},
    {"_AlphaSimR_callRRBLUP_GCA2",   (DL_FUNC) &_AlphaSimR_callRRBLUP_GCA2,   13},
    {"_AlphaSimR_callRRBLUP_MV",     (DL_FUNC) &_AlphaSimR_callRRBLUP_MV,      8},
    {"_AlphaSimR_callRRBLUP_SCA",    (DL_FUNC) &_AlphaSimR_callRRBLUP_SCA,     8},
    {"_AlphaSimR_callRRBLUP2",       (DL_FUNC) &_AlphaSimR_callRRBLUP2,       12},
    {"_AlphaSimR_createDH2",         (DL_FUNC) &_AlphaSimR_createDH2,          5},
    {"_AlphaSimR_cross",             (DL_FUNC) &_AlphaSimR_cross,             13},
    {"_AlphaSimR_getGeno",           (DL_FUNC) &_AlphaSimR_getGeno,            4},
    {"_AlphaSimR_getGv",             (DL_FUNC) &_AlphaSimR_getGv,              3},
    {"_AlphaSimR_getHaplo",          (DL_FUNC) &_AlphaSimR_getHaplo,           4},
    {"_AlphaSimR_getHybridGv",       (DL_FUNC) &_AlphaSimR_getHybridGv,        6},
    {"_AlphaSimR_getIbdHaplo",       (DL_FUNC) &_AlphaSimR_getIbdHaplo,        3},
    {"_AlphaSimR_getIbdRecHist",     (DL_FUNC) &_AlphaSimR_getIbdRecHist,      3},
    {"_AlphaSimR_getMaternalGeno",   (DL_FUNC) &_AlphaSimR_getMaternalGeno,    4},
    {"_AlphaSimR_getNumThreads",     (DL_FUNC) &_AlphaSimR_getNumThreads,      0},
    {"_AlphaSimR_getOneHaplo",       (DL_FUNC) &_AlphaSimR_getOneHaplo,        5},
    {"_AlphaSimR_getPaternalGeno",   (DL_FUNC) &_AlphaSimR_getPaternalGeno,    4},
    {"_AlphaSimR_MaCS",              (DL_FUNC) &_AlphaSimR_MaCS,               5},
    {"_AlphaSimR_mergeGeno",         (DL_FUNC) &_AlphaSimR_mergeGeno,          2},
    {"_AlphaSimR_mergeMultGeno",     (DL_FUNC) &_AlphaSimR_mergeMultGeno,      4},
    {"_AlphaSimR_mergeMultIntMat",   (DL_FUNC) &_AlphaSimR_mergeMultIntMat,    3},
    {"_AlphaSimR_packHaplo",         (DL_FUNC) &_AlphaSimR_packHaplo,          3},
    {"_AlphaSimR_popVar",            (DL_FUNC) &_AlphaSimR_popVar,             1},
    {"_AlphaSimR_sampAllComb",       (DL_FUNC) &_AlphaSimR_sampAllComb,        3},
    {"_AlphaSimR_sampHalfDialComb",  (DL_FUNC) &_AlphaSimR_sampHalfDialComb,   2},
    {"_AlphaSimR_writeASGenotypes",  (DL_FUNC) &_AlphaSimR_writeASGenotypes,   7},
    {"_AlphaSimR_writeASHaplotypes", (DL_FUNC) &_AlphaSimR_writeASHaplotypes,  7},
    {"_AlphaSimR_writeGeno",         (DL_FUNC) &_AlphaSimR_writeGeno,          5},
    {"_AlphaSimR_writeOneHaplo",     (DL_FUNC) &_AlphaSimR_writeOneHaplo,      6},
    {"_AlphaSimR_writePlinkPed",     (DL_FUNC) &_AlphaSimR_writePlinkPed,      5},
    {NULL, NULL, 0}
};

void R_init_AlphaSimR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
