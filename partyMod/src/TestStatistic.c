
/**
    Test statistics for conditional inference
    *\file TestStatistic.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"


/**
    Standardizes a statistic t of length pq with mean mu and covariance Sigma
    for variances > tol \n
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_standardize(const double *t, const double *mu, const double *Sigma, 
                   int pq, double tol, double *ans) {
                   
    int i;
    double sd;
    
    for (i = 0; i < pq; i++) {
        sd = Sigma[i*pq + i]; 
        if (sd > tol)
            ans[i] = (t[i] - mu[i])/sqrt(sd);
        else
            ans[i] = 0.0;
    }
}


/**
    Absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_absstandardize(const double *t, const double *mu, const double *Sigma,
                      int pq, double tol, double *ans) {
                      
    C_standardize(t, mu, Sigma, pq, tol, ans);
    C_abs(ans, pq);
}


/**
    Maximum absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
*/

double C_maxabsTestStatistic(const double *t, const double *mu, const double *Sigma,
                             int pq, double tol) {
                           
     double *mem, ans;
     
     mem = Calloc(pq, double);
     C_absstandardize(t, mu, Sigma, pq, tol, mem);
     ans = C_max(mem, pq);
     Free(mem);
     return(ans);
}


/**
    R-interface to C_maxabsTestStatistic
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param tol tolerance for variances
*/

SEXP R_maxabsTestStatistic(SEXP t, SEXP mu, SEXP Sigma, SEXP tol) {

     SEXP ans;
     int pq;
     
     pq = LENGTH(t);
     
     PROTECT(ans = allocVector(REALSXP, 1));
     REAL(ans)[0] = C_maxabsTestStatistic(REAL(t), REAL(mu), REAL(Sigma), pq, 
                                          REAL(tol)[0]);
     UNPROTECT(1);
     return(ans);
}


/**
    Quadratic form t(t - mu) SigmaPlus (t - mu) \n
    *\param t the vector of statistics
    *\param mu expectations
    *\param SigmaPlus Moore-Penrose inverse
    *\param pq dimension of t
*/
                                
double C_quadformTestStatistic(const double *t, const double *mu, 
                               const double *SigmaPlus, int pq) {

    int i, j;
    double quadform = 0.0, *tmmu, *tmmuSigmaPlus;

    tmmu = Calloc(pq, double);
    for (i = 0; i < pq; i++)
        tmmu[i] = t[i] - mu[i];
    
    tmmuSigmaPlus = Calloc(pq, double);
    for (i = 0; i < pq; i++)  {
        tmmuSigmaPlus[i] = 0.0;
        for (j = 0; j < pq; j++)
            tmmuSigmaPlus[i] += tmmu[j] * SigmaPlus[i * pq + j];
        quadform += tmmuSigmaPlus[i] * tmmu[i];
    }

    Free(tmmu); Free(tmmuSigmaPlus);
    return(quadform);
}


/**
    R-interface to C_quadformTestStatistic \n
    *\param t the vector of statistics
    *\param mu expectations
    *\param SigmaPlus Moore-Penrose inverse
*/
    
SEXP R_quadformTestStatistic(SEXP t, SEXP mu, SEXP SigmaPlus) {

    SEXP ans;
    int pq;
    
    pq = LENGTH(t);
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = C_quadformTestStatistic(REAL(t), 
                                           REAL(mu), REAL(SigmaPlus), pq);
    UNPROTECT(1);
    return(ans);
}
