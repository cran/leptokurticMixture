// See section 6.6.1 of Writing R Extensions (Fortran character strings)
#define USE_FC_LEN_T

#include <R_ext/Lapack.h>
#include <stdbool.h>

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP riccatiCareSolution(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"riccatiCareSolution", (DL_FUNC) &riccatiCareSolution, 1},
    {NULL, NULL, 0}
};

void R_init_leptokurticMixture(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



#ifndef FCONE
# define FCONE
#endif
int sortbynegreal(const double* real, const double* img) {
  if (*real < 0  )
    return(true);
  else
    return(false);
}

SEXP riccatiCareSolution(SEXP x) {
  int *dim, n, m, info, izero = 0, lwork = -1;
  double *work, tmp;
  
  dim = INTEGER(getAttrib(x, R_DimSymbol));
  
  n = dim[0];
  int bwork[n];
  m = n/2;
  
  dgees_("V", "S", sortbynegreal, dim, (double *) NULL, dim, &izero,
         (double *) NULL, (double *) NULL, (double *) NULL, dim,
         &tmp, &lwork, bwork, &info FCONE FCONE);
  
  lwork = (int) tmp;
  work = (double*)malloc( lwork*sizeof(double) );
  
  double *wr, *wi, *matz;
  wr = (double*)malloc( n*sizeof(double) );
  wi = (double*)malloc( n*sizeof(double) );
  matz = (double*)malloc( n*n*sizeof(double) );
  
  dgees_("V", "S", sortbynegreal, dim,  REAL(x), dim, &izero, wr, wi,
         matz, dim, work, &lwork, bwork, &info FCONE FCONE);
  
  SEXP ans = PROTECT(allocMatrix(REALSXP, m, m));
  
  double *px = REAL(ans), *py;
  py = (double*)malloc( m*m*sizeof(double) );
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      px[i*m+j] = matz[ i + j*n + m ];
      py[i*m+j] = matz[ i + j*n ];
    }
  }
  
  int *ipiv;
  ipiv = (int*)malloc( m*sizeof(int) );
  int num = m, nrhs = m, lda = m, ldb = m, info2;
  
  dgesv_( &num, &nrhs, py, &lda, ipiv, px, &ldb, &info2);
  
  free(matz); free(wr); free(wi); free(ipiv);  free(work);
  
  UNPROTECT(1);
  return ans;
}


