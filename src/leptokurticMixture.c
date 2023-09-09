

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#ifndef ILP64
# define La_INT int
# define La_LGL int
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdbool.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>


extern SEXP riccatiCareSolution(SEXP x);

static const R_CallMethodDef CallEntries[] = {
    {"riccatiCareSolution", (DL_FUNC) &riccatiCareSolution, 1},
    {NULL, NULL, 0}
};

void R_init_leptokurticMixture(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


//La_extern La_LGL F77_NAME(selector)(const double* real, const double* img );

La_LGL sortbynegreal(const double* real, const double* img) {
  if (*real < 0  )
    return(true);
  else
    return(false);
}



extern SEXP riccatiCareSolution(SEXP x) {
  int *dim, n, m, info, izero = 0, lwork = -1;
  double *work, tmp;
  
  dim = INTEGER(getAttrib(x, R_DimSymbol));

  n = dim[0];
  int bwork[n];
  m = n/2;
  
  F77_CALL(dgees)("V", "S", &sortbynegreal, dim, (double *) NULL, dim, &izero,
            (double *) NULL, (double *) NULL, (double *) NULL, dim,
            &tmp, &lwork, bwork, &info FCONE FCONE);
    
  lwork = (int) tmp;
  work = (double*)malloc( lwork*sizeof(double) );
  
  double *wr, *wi, *matz;
  wr = (double*)malloc( n*sizeof(double) );
  wi = (double*)malloc( n*sizeof(double) );
  matz = (double*)malloc( n*n*sizeof(double) );

  F77_CALL(dgees)("V", "S", &sortbynegreal, dim,  REAL(x), dim, &izero, wr, wi,
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
  
  F77_CALL(dgesv)( &num, &nrhs, py, &lda, ipiv, px, &ldb, &info2 );
    
  free(matz); free(wr); free(wi); free(ipiv);  free(work);
  
  UNPROTECT(1);
  return ans;
}




