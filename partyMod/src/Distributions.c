
/**
    Conditional Distributions
    *\file Distributions.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"


/**
    Conditional asymptotic P-value of a quadratic form\n
    *\param tstat test statistic
    *\param df degree of freedom
*/
        
double C_quadformConditionalPvalue(const double tstat, const double df) {
    return(pchisq(tstat, df, 0, 0));
}


/**
    R-interface to C_quadformConditionalPvalue\n
    *\param tstat test statitstic
    *\param df degree of freedom
*/

SEXP R_quadformConditionalPvalue(SEXP tstat, SEXP df) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = C_quadformConditionalPvalue(REAL(tstat)[0], REAL(df)[0]);
    UNPROTECT(1);
    return(ans);
}


/**
    Conditional asymptotic P-value of a maxabs-type test statistic\n
    Basically the functionality from package `mvtnorm' \n
    *\param tstat test statitstic
    *\param Sigma covariance matrix
    *\param pq nrow(Sigma)
    *\param maxpts number of Monte-Carlo steps
    *\param releps relative error
    *\param abseps absolute error
    *\param tol tolerance
*/

double C_maxabsConditionalPvalue(const double tstat, const double *Sigma, 
    const int pq, int *maxpts, double *releps, double *abseps, double *tol) {

    int *n, *nu, *inform, i, j, *infin, sub, *index, nonzero, iz, jz;
    double *lower, *upper, *delta, *corr, *sd, *myerror,
           *prob, ans;

    /* univariate problem */
    if (pq == 1) 
        return(2*pnorm(fabs(tstat)*-1.0, 0.0, 1.0, 1, 0)); /* return P-value */
    
    n = Calloc(1, int);
    nu = Calloc(1, int);
    myerror = Calloc(1, double);
    prob = Calloc(1, double);
    nu[0] = 0;
    inform = Calloc(1, int);
    n[0] = pq;
        
    if (n[0] == 2)  
         corr = Calloc(1, double);
    else 
         corr = Calloc(n[0] + ((n[0] - 2) * (n[0] - 1))/2, double);
    
    sd = Calloc(n[0], double);
    lower = Calloc(n[0], double);
    upper = Calloc(n[0], double);
    infin = Calloc(n[0], int);
    delta = Calloc(n[0], double);
    index = Calloc(n[0], int);

    /* determine elements with non-zero variance */ 

    nonzero = 0;
    for (i = 0; i < n[0]; i++) {
        if (Sigma[i*n[0] + i] > tol[0]) {
            index[nonzero] = i;
            nonzero++;
        }
    }

    /* mvtdst assumes the unique elements of the triangular 
       covariance matrix to be passes as argument CORREL 
    */

    for (iz = 0; iz < nonzero; iz++) {

        /* handle elements with non-zero variance only */
        i = index[iz];

        /* standard deviations */
        sd[i] = sqrt(Sigma[i*n[0] + i]);
                
        /* always look at the two-sided problem */           
        lower[iz] = fabs(tstat) * -1.0;
        upper[iz] = fabs(tstat);
        infin[iz] = 2;
        delta[iz] = 0.0;
        
        /* set up vector of correlations, i.e., the upper 
           triangular part of the covariance matrix) */
        for (jz = 0; jz < iz; jz++) {
            j = index[jz]; 
            sub = (int) (jz + 1) + (double) ((iz - 1) * iz) / 2 - 1;
            if (sd[i] == 0.0 || sd[j] == 0.0) 
                corr[sub] = 0.0; 
            else 
                corr[sub] = Sigma[i*n[0] + j] / (sd[i] * sd[j]);
        }
    }
    n[0] = nonzero;
        
    /* call FORTRAN subroutine */
    F77_CALL(mvtdst)(n, nu, lower, upper, infin, corr, delta, 
                     maxpts, abseps, releps, myerror, prob, inform);
                         
    /* inform == 0 means: everything is OK */
    switch (inform[0]) {
        case 0: break;
        case 1: warning("cmvnorm: completion with ERROR > EPS"); break;
        case 2: warning("cmvnorm: N > 1000 or N < 1"); 
                prob[0] = 0.0; 
                break;
        case 3: warning("cmvnorm: correlation matrix not positive semi-definite"); 
                prob[0] = 0.0; 
                break;
        default: warning("cmvnorm: unknown problem in MVTDST");
                 prob[0] = 0.0;
    }
    ans = prob[0];
    Free(corr); Free(sd); Free(lower); Free(upper); 
    Free(infin); Free(delta); Free(myerror); Free(prob);
    Free(n); Free(nu); Free(inform); 
    return(1 - ans);  /* return P-value */
}


/**
   R-interface to C_maxabsConditionalPvalue \n
   *\param tstat test statitstic
   *\param Sigma covariance matrix
   *\param maxpts number of Monte-Carlo steps
   *\param releps relative error
   *\param abseps absolute error
   *\param tol tolerance
*/

SEXP R_maxabsConditionalPvalue(SEXP tstat, SEXP Sigma, SEXP maxpts, 
                               SEXP releps, SEXP abseps, SEXP tol) {
           
    SEXP ans;                    
    int pq;
    
    pq = nrow(Sigma);
   
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = C_maxabsConditionalPvalue(REAL(tstat)[0], REAL(Sigma), pq, 
        INTEGER(maxpts), REAL(releps), REAL(abseps), REAL(tol));
    UNPROTECT(1);
    return(ans);
}


/**
    Monte-Carlo approximation to the conditional pvalues 
    *\param criterion vector of node criteria for each input 
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param varctrl an object of class `VariableControl'
    *\param gtctrl an object of class `GlobalTestControl'
    *\param ans_pvalues return values; vector of adjusted pvalues 
*/

void C_MonteCarlo(double *criterion, SEXP learnsample, SEXP weights, 
                  SEXP fitmem, SEXP varctrl, SEXP gtctrl, double *ans_pvalues) {

    int ninputs, nobs, j, i, k;
    SEXP responses, inputs, y, x, xmem, expcovinf;
    double sweights, *stats, tmp = 0.0, smax, *dweights; 
    int m, *counts, b, B, *dummy, *permindex, *index, *permute;
    
    ninputs = get_ninputs(learnsample);
    nobs = get_nobs(learnsample);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    dweights = REAL(weights);
    
    /* number of Monte-Carlo replications */
    B = get_nresample(gtctrl);
    
    /* y = get_transformation(responses, 1); */
    y = get_test_trafo(responses);
    
    expcovinf = GET_SLOT(fitmem, PL2_expcovinfSym);

    sweights = REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0];
    m = (int) sweights;
    
    stats = Calloc(ninputs, double);
    counts = Calloc(ninputs, int);
    
    dummy = Calloc(m, int);
    permute = Calloc(m, int);
    index = Calloc(m, int);
    permindex = Calloc(m, int);
                
    /* expand weights, see appendix of 
       `Unbiased Recursive Partitioning: A Conditional Inference Framework' */
    j = 0;
    for (i = 0; i < nobs; i++) {
        for (k = 0; k < dweights[i]; k++) {
            index[j] = i;
            j++;
        }
    }

    for (b = 0; b < B; b++) {

        /* generate a admissible permutation */
        C_SampleNoReplace(dummy, m, m, permute);
        for (k = 0; k < m; k++) permindex[k] = index[permute[k]];

        /* for all input variables */
        for (j = 1; j <= ninputs; j++) {
            x = get_transformation(inputs, j);

            /* compute test statistic or pvalue for the permuted data */
            xmem = get_varmemory(fitmem, j);
            if (!has_missings(inputs, j)) {
                C_PermutedLinearStatistic(REAL(x), ncol(x), REAL(y), ncol(y), 
                    nobs, m, index, permindex, 
                    REAL(GET_SLOT(xmem, PL2_linearstatisticSym)));
            } else {
                error("cannot resample with missing values");
            }
            
            /* compute the criterion, i.e. something to be MAXIMISED */
            C_TeststatCriterion(xmem, varctrl, &tmp, &stats[j - 1]);
        }
        
        /* the maximum of the permuted test statistics / 1 - pvalues */
        smax = C_max(stats, ninputs);

        /* count the number of permuted > observed */
        for (j = 0; j < ninputs; j++) {
            if (smax > criterion[j]) counts[j]++;
        }
    }
    
    /* return adjusted pvalues */
    for (j = 0; j < ninputs; j++)
        ans_pvalues[j] = (double) counts[j] / B;
                
    /* <FIXME> we try to assess the linear statistics later on 
               (in C_Node, for categorical variables) 
               but have used this memory for resampling here */

    for (j = 1; j <= ninputs; j++) {
        x = get_transformation(inputs, j);
        /* re-compute linear statistics for unpermuted data */
        xmem = get_varmemory(fitmem, j);
        C_LinearStatistic(REAL(x), ncol(x), REAL(y), ncol(y), 
                      dweights, nobs, 
                      REAL(GET_SLOT(xmem, PL2_linearstatisticSym)));
    }
    /* </FIXME> */
    
    Free(stats); Free(counts); Free(dummy); Free(permute); 
    Free(index); Free(permindex);
}


/**
    R-interface to C_MonteCarlo \n
    *\param criterion vector of node criteria for each input
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param varctrl an object of class `VariableControl'
    *\param gtctrl an object of class `GlobalTestControl'
*/

SEXP R_MonteCarlo(SEXP criterion, SEXP learnsample, SEXP weights, 
                  SEXP fitmem, SEXP varctrl, SEXP gtctrl) {
                  
     SEXP ans;
     
     GetRNGstate();
     
     PROTECT(ans = allocVector(REALSXP, get_ninputs(learnsample)));
     C_MonteCarlo(REAL(criterion), learnsample, weights, fitmem, varctrl, 
                  gtctrl, REAL(ans));
                  
     PutRNGstate();
                  
     UNPROTECT(1);
     return(ans);
}
