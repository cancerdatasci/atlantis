
double C_quadformConditionalPvalue(const double tstat, const double df);
double C_maxabsConditionalPvalue(const double tstat, const double *Sigma, const int pq,
                                 int *maxpts, double *releps, double *abseps, double *tol);
void C_MonteCarlo(double *pvalues, SEXP learnsample, SEXP weights,
                  SEXP fitmem, SEXP varctrl, SEXP gtctrl, double *ans);
