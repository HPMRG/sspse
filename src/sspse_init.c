#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .C calls */
extern void gcmp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gnbinom(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gpln(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void lcmp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rcmp(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
	    {"gcmp",            (DL_FUNC) &gcmp,            23},
	    {"gnbinom",         (DL_FUNC) &gnbinom,         21},
	    {"gpln",            (DL_FUNC) &gpln,            21},
	    {"lcmp",            (DL_FUNC) &lcmp,            23},
	    {"rcmp",            (DL_FUNC) &rcmp,             6},
	    {NULL, NULL, 0}
};

void R_init_sspse(DllInfo *dll)
{
	    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	        R_useDynamicSymbols(dll, FALSE);
}
