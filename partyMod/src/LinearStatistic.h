
void C_LinearStatistic(const double *x, const int p,
                       const double *y, const int q,
                       const double *weights, const int n,
                       double *ans);
void C_ExpectCovarInfluence(const double* y, const int q,
                            const double* weights, const int n,
                            SEXP ans);
void C_ExpectCovarLinearStatistic(const double* x, const int p,
                                  const double* y, const int q,
                                  const double* weights, const int n,
                                  const SEXP expcovinf, SEXP ans);
void C_PermutedLinearStatistic(const double *x, const int p,
                               const double *y, const int q,
                               const int n, const int nperm,
                               const int *indx, const int *perm,
                               double *ans);
SEXP R_ExpectCovarInfluence(SEXP y, SEXP weights);
