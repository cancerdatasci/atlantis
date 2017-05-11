
void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans);
/* SEXP La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, 
            SEXP v, SEXP method); */
void C_SampleNoReplace(int *x, int m, int k, int *ans);
void C_MPinv (SEXP x, double tol, SEXP svdmem, SEXP ans);
double C_max(const double *x, const int n);
void C_abs(double *x, int n);
void C_matprod(double *x, int nrx, int ncx,
               double *y, int nry, int ncy, double *z);
void C_matprodT(double *x, int nrx, int ncx,
                double *y, int nry, int ncy, double *z);
int nrow(SEXP x);
int ncol(SEXP y);
int C_whichmax(double *pvalue, double *teststat, int ninputs);
int i_in_set(int i, int *iset, int p);
int C_i_in_set(int i, SEXP set);
void C_SampleSplitting(int n, double *prob, int *weights, int k);
void C_remove_weights(SEXP subtree, int removestats);
double* C_tempweights(int j, SEXP weights, SEXP fitmem, SEXP inputs);
void C_linexpcovReduce (SEXP x);
