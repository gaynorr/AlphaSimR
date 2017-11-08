#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP AlphaSimR_calcChrFreq(SEXP);
extern SEXP AlphaSimR_calcG(SEXP);
extern SEXP AlphaSimR_calcGenParam(SEXP, SEXP);
extern SEXP AlphaSimR_calcGIbs(SEXP);
extern SEXP AlphaSimR_callRRBLUP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_callRRBLUP_GCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_callRRBLUP_MV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_callRRBLUP_SCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_changeId(SEXP, SEXP);
extern SEXP AlphaSimR_convToImat(SEXP);
extern SEXP AlphaSimR_createDH2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_cross2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_crossPedigree(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_fastDist(SEXP);
extern SEXP AlphaSimR_fastPairDist(SEXP, SEXP);
extern SEXP AlphaSimR_gaussKernel(SEXP, SEXP);
extern SEXP AlphaSimR_gebvGCA(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_gebvRR(SEXP, SEXP);
extern SEXP AlphaSimR_gebvSCA(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_getDomGeno(SEXP);
extern SEXP AlphaSimR_getGeno(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_getGv(SEXP, SEXP);
extern SEXP AlphaSimR_getHaplo(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_getHybridGv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_getOneHaplo(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_MaCS(SEXP, SEXP);
extern SEXP AlphaSimR_mergeGeno(SEXP, SEXP);
extern SEXP AlphaSimR_packHaplo(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_popVar(SEXP);
extern SEXP AlphaSimR_readMat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_sampAllComb(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_sampHalfDialComb(SEXP, SEXP);
extern SEXP AlphaSimR_solveMKM(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_solveMVM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_solveRRBLUP(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_solveRRBLUPMK(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_solveRRBLUPMV(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_solveUVM(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_tuneTraitA(SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_tuneTraitAD(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_writeASGenotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_writeASHaplotypes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_writeGeno(SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_writeOneHaplo(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AlphaSimR_zero();

static const R_CallMethodDef CallEntries[] = {
    {"AlphaSimR_calcChrFreq",       (DL_FUNC) &AlphaSimR_calcChrFreq,       1},
    {"AlphaSimR_calcG",             (DL_FUNC) &AlphaSimR_calcG,             1},
    {"AlphaSimR_calcGenParam",      (DL_FUNC) &AlphaSimR_calcGenParam,      2},
    {"AlphaSimR_calcGIbs",          (DL_FUNC) &AlphaSimR_calcGIbs,          1},
    {"AlphaSimR_callRRBLUP",        (DL_FUNC) &AlphaSimR_callRRBLUP,        6},
    {"AlphaSimR_callRRBLUP_GCA",    (DL_FUNC) &AlphaSimR_callRRBLUP_GCA,    7},
    {"AlphaSimR_callRRBLUP_MV",     (DL_FUNC) &AlphaSimR_callRRBLUP_MV,     7},
    {"AlphaSimR_callRRBLUP_SCA",    (DL_FUNC) &AlphaSimR_callRRBLUP_SCA,    7},
    {"AlphaSimR_changeId",          (DL_FUNC) &AlphaSimR_changeId,          2},
    {"AlphaSimR_convToImat",        (DL_FUNC) &AlphaSimR_convToImat,        1},
    {"AlphaSimR_createDH2",         (DL_FUNC) &AlphaSimR_createDH2,         5},
    {"AlphaSimR_cross2",            (DL_FUNC) &AlphaSimR_cross2,            6},
    {"AlphaSimR_crossPedigree",     (DL_FUNC) &AlphaSimR_crossPedigree,     4},
    {"AlphaSimR_fastDist",          (DL_FUNC) &AlphaSimR_fastDist,          1},
    {"AlphaSimR_fastPairDist",      (DL_FUNC) &AlphaSimR_fastPairDist,      2},
    {"AlphaSimR_gaussKernel",       (DL_FUNC) &AlphaSimR_gaussKernel,       2},
    {"AlphaSimR_gebvGCA",           (DL_FUNC) &AlphaSimR_gebvGCA,           4},
    {"AlphaSimR_gebvRR",            (DL_FUNC) &AlphaSimR_gebvRR,            2},
    {"AlphaSimR_gebvSCA",           (DL_FUNC) &AlphaSimR_gebvSCA,           3},
    {"AlphaSimR_getDomGeno",        (DL_FUNC) &AlphaSimR_getDomGeno,        1},
    {"AlphaSimR_getGeno",           (DL_FUNC) &AlphaSimR_getGeno,           3},
    {"AlphaSimR_getGv",             (DL_FUNC) &AlphaSimR_getGv,             2},
    {"AlphaSimR_getHaplo",          (DL_FUNC) &AlphaSimR_getHaplo,          3},
    {"AlphaSimR_getHybridGv",       (DL_FUNC) &AlphaSimR_getHybridGv,       5},
    {"AlphaSimR_getOneHaplo",       (DL_FUNC) &AlphaSimR_getOneHaplo,       4},
    {"AlphaSimR_MaCS",              (DL_FUNC) &AlphaSimR_MaCS,              2},
    {"AlphaSimR_mergeGeno",         (DL_FUNC) &AlphaSimR_mergeGeno,         2},
    {"AlphaSimR_packHaplo",         (DL_FUNC) &AlphaSimR_packHaplo,         3},
    {"AlphaSimR_popVar",            (DL_FUNC) &AlphaSimR_popVar,            1},
    {"AlphaSimR_readMat",           (DL_FUNC) &AlphaSimR_readMat,           6},
    {"AlphaSimR_sampAllComb",       (DL_FUNC) &AlphaSimR_sampAllComb,       3},
    {"AlphaSimR_sampHalfDialComb",  (DL_FUNC) &AlphaSimR_sampHalfDialComb,  2},
    {"AlphaSimR_solveMKM",          (DL_FUNC) &AlphaSimR_solveMKM,          4},
    {"AlphaSimR_solveMVM",          (DL_FUNC) &AlphaSimR_solveMVM,          6},
    {"AlphaSimR_solveRRBLUP",       (DL_FUNC) &AlphaSimR_solveRRBLUP,       3},
    {"AlphaSimR_solveRRBLUPMK",     (DL_FUNC) &AlphaSimR_solveRRBLUPMK,     3},
    {"AlphaSimR_solveRRBLUPMV",     (DL_FUNC) &AlphaSimR_solveRRBLUPMV,     5},
    {"AlphaSimR_solveUVM",          (DL_FUNC) &AlphaSimR_solveUVM,          4},
    {"AlphaSimR_tuneTraitA",        (DL_FUNC) &AlphaSimR_tuneTraitA,        3},
    {"AlphaSimR_tuneTraitAD",       (DL_FUNC) &AlphaSimR_tuneTraitAD,       5},
    {"AlphaSimR_writeASGenotypes",  (DL_FUNC) &AlphaSimR_writeASGenotypes,  7},
    {"AlphaSimR_writeASHaplotypes", (DL_FUNC) &AlphaSimR_writeASHaplotypes, 7},
    {"AlphaSimR_writeGeno",         (DL_FUNC) &AlphaSimR_writeGeno,         4},
    {"AlphaSimR_writeOneHaplo",     (DL_FUNC) &AlphaSimR_writeOneHaplo,     5},
    {"AlphaSimR_zero",              (DL_FUNC) &AlphaSimR_zero,              0},
    {NULL, NULL, 0}
};

void R_init_AlphaSimR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
