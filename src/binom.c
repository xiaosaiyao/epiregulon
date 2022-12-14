#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <R_ext/Rdynload.h>


SEXP to_start(SEXP start_val, SEXP n, SEXP start_k, SEXP prob){
  int k_st, N, k;
  double max_prob, p_quot;
  double *vec_p;
  SEXP vec;
  N = INTEGER(n)[0];
  k_st = INTEGER(start_k)[0];
  max_prob = REAL(start_val)[0];
  p_quot = (1-REAL(prob)[0])/REAL(prob)[0];
  k = k_st;
  PROTECT(vec = Rf_allocVector(REALSXP, k_st + 1));
  vec_p = REAL(vec);
  vec_p[0] = max_prob;
  for(int i = 1; i <= k_st; ++i){
    vec_p[i] = vec_p[i-1]*(k)/(N-k+1)*p_quot;
    k-=1;
  }
  UNPROTECT(1);
  return(vec);
  }


SEXP to_end(SEXP start_val, SEXP n, SEXP start_k, SEXP prob){
    int k_st, N, k;
    double max_prob, p_quot;
    double *vec_p;
    SEXP vec;
    N = INTEGER(n)[0];
    k_st = INTEGER(start_k)[0];
    max_prob = REAL(start_val)[0];
    k = k_st;
    p_quot = REAL(prob)[0]/(1-REAL(prob)[0]);
    PROTECT(vec = Rf_allocVector(REALSXP, N + 2 - k_st));
    vec_p = REAL(vec);
    vec_p[0] = max_prob;
    for(int i = 1; i<=(N + 1 - k_st); i++){
      vec_p[i] = vec_p[i-1]*(N-k)/(k+1)*p_quot;
      k+=1;
    }
    UNPROTECT(1);
    return(vec);
  }

static const R_CallMethodDef callMethods[]  = {
  {"to_start", (DL_FUNC) &to_start, 4},
  {"to_end", (DL_FUNC) &to_end, 4},
  {NULL, NULL, 0}
};

void
  R_init_binom(DllInfo *info)
  {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
  }
