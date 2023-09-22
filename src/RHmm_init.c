#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP RBaumWelch(SEXP, SEXP, SEXP);
extern SEXP RComputeCov(SEXP, SEXP);
extern SEXP Rforwardbackward(SEXP, SEXP, SEXP);
extern SEXP RScoreAndInformation(SEXP, SEXP);
extern SEXP RViterbi(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"RBaumWelch",           (DL_FUNC) &RBaumWelch,           3},
    {"RComputeCov",          (DL_FUNC) &RComputeCov,          2},
    {"Rforwardbackward",     (DL_FUNC) &Rforwardbackward,     3},
    {"RScoreAndInformation", (DL_FUNC) &RScoreAndInformation, 2},
    {"RViterbi",             (DL_FUNC) &RViterbi,             2},
    {NULL, NULL, 0}
};

void R_init_RHmm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}