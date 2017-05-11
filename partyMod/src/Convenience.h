                
void C_LinStatExpCov(const double *x, const int p,
                     const double *y, const int q,
                     const double *weights, const int n,
                     const int cexpcovinf, SEXP expcovinf, SEXP ans);
void C_LinStatExpCovMPinv(SEXP linexpcov, double tol);
void C_MLinearStatistic(SEXP linexpcov, SEXP ScoreMatrix, SEXP ans);
double C_TestStatistic(const SEXP linexpcov, const int type, const double tol);
double C_ConditionalPvalue(const double tstat, SEXP linexpcov,
                           const int type, double tol,
                           int *maxpts, double *releps, double *abseps);
