
void C_standardize(const double *t, const double *mu, const double *Sigma,
                   int pq, double tol, double *ans);
double C_maxabsTestStatistic(const double *t, const double *mu, const double *Sigma,
                             int pq, double tol);
double C_quadformTestStatistic(const double *t, const double *mu,
                               const double *SigmaPlus, int pq);
