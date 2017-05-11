
/**
    Binary splits 
    *\file Splits.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"


/**
    Search for a cutpoint in a ordered variable x maximizing a two-sample
    statistic w.r.t. (the influence function of ) the response variable y.
    *\param x raw numeric measurements 
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param orderx the ordering of the transformations, i.e. R> order(x)
    *\param splitctrl an object of class `SplitControl'
    *\param linexpcov2sample an (uninitialized) object of class 
                             `LinStatExpectCovar' with p = 1
    *\param expcovinf an initialized object of class `ExpectCovarInfluence'
    *\param cutpoint return value; pointer to a double for the cutpoint in x
    *\param maxstat return value; pointer to a double for the 
                    maximal test statistic
    *\param statistics return value; pointer to a n-dim double 
                       for the statistics
*/
                                
void C_split(const double *x, int p,
             const double *y, int q,
             const double *weights, int n,
             const int *orderx,
             SEXP splitctrl, SEXP linexpcov2sample, 
             SEXP expcovinf, double *cutpoint, double *maxstat, 
             double *statistics) {

    double *dExp_y, *dCov_y, *dlinstat, *dexpect, *dcovar, 
           tol, sweights, minprob, minbucket, w, tx, f1, f2, f1w, f2ww, tmp;
    double minobs, maxobs, xmax;
    int lastj, i, j, k;

    if (p != 1) error("C_split: p not equal to one");
    tol = get_tol(splitctrl);

    /* init statistics and determine the maximal value with positive weight 
       since we can't choose this one as cutpoint
    */
    xmax = 0.0;
    for (i = 0; i < n; i++) {
        statistics[i] = 0.0;
        if (weights[i] > 0.0 && x[i] > xmax) xmax = x[i];
    }

    /* we already have expecation and covariance of the response
     * values and the sum of the weights */
    dExp_y = REAL(GET_SLOT(expcovinf, PL2_expectationSym));
    dCov_y = REAL(GET_SLOT(expcovinf, PL2_covarianceSym));
    sweights = REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0];

    /* if there is something to split */
    if (sweights > 1) {

        /* we need to ensure that at least minbucket weights 
           are there to split (either left or right) */
        minprob = get_minprob(splitctrl);
        minbucket = get_minbucket(splitctrl);
        minobs = sweights * minprob + 1.0;

        if (minobs < minbucket) 
            minobs = minbucket; 
        maxobs = sweights * (1 - minprob) - 1.0;
        if (maxobs > sweights - minbucket) 
            maxobs = sweights - minbucket; 

        f1 = (double) sweights / (sweights - 1);
        f2 = 1.0 / (sweights - 1);
        w = 0.0;
    
        /* pointers to the R-objects */
        dlinstat = REAL(GET_SLOT(linexpcov2sample, PL2_linearstatisticSym));
        for (k = 0; k < q; k++) dlinstat[k] = 0.0;
        dexpect = REAL(GET_SLOT(linexpcov2sample, PL2_expectationSym));
        dcovar = REAL(GET_SLOT(linexpcov2sample, PL2_covarianceSym));

        tx = 0.0;
        lastj = 0;

        /* for all possible cutpoints (defined by the observations x) */
        for (i = 0; i < (n - 1); i++) {
    
            /* the ordering of the ith observation */
            j = orderx[i] - 1;
        
            /* if the corresponding weight is zero */
            if (weights[j] == 0.0) continue;

            /* just a check: can be removed later */
            if (w > 0 && x[j] < tx)
                warning("C_split: inconsistent ordering: %f < %f!\n", 
                        x[j], tx);
        
            /* handle ties: delete the entry of the last visited observation
               (take care of zero weights!) */
            if (w > 0 && x[j] == tx)
                statistics[lastj] = 0.0; 

            /* store the value and position of the j smallest observation */
            tx = x[j];
            lastj = j;
        
            w += weights[j];

            /* do not consider those splits */
            if (w > maxobs || x[j] >= xmax) break;

            /* compute the linear statistic and expectation and 
             * covariance if needed */
            for (k = 0; k < q; k++)
                dlinstat[k] += y[n * k + j] * weights[j];
 
            if (w >= minobs) {
                for (k = 0; k < q; k++)
                    dexpect[k] = w * dExp_y[k];

                f1w = f1 * w;
                f2ww = f2 * w * w;
                for (k = 0; k < q*q; k++)
                    dcovar[k] = f1w * dCov_y[k] - f2ww * dCov_y[k];
            } else {
                continue;
            }
        
            /* the absolute standardized test statistic, to be maximized */
            /* statistics[j] = C_maxabsTestStatistic(dlinstat, 
                   dexpect, dcovar, q, tol); */

            /* much faster but uses maxabs always*/
            statistics[j] = 0.0;
            for (k = 0; k < q; k++) {
                if (dcovar[k * q + k] <= tol) continue;
                tmp = fabs(dlinstat[k] - dexpect[k]) / sqrt(dcovar[k * q + k]);
                if (statistics[j] < tmp) statistics[j] = tmp;
            }

        }
    
        /* search for the maximum and the best separating cutpoint */
        /* <FIXME> the result might differ between 32 and 64bit systems 
                   because of rounding errors in 'statistics' */
        maxstat[0] = 0.0;        
        for (i = 0; i < n; i++) {
            if (statistics[i] > maxstat[0]) {
                maxstat[0] = statistics[i];
                cutpoint[0] = x[i];
            }
        }
        /* </FIXME> */
    }
}


/**
    R-interface to C_split (does not handle ordered y's)
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param orderx the ordering of the transformations, i.e. R> order(x)
    *\param linexpcov2sample an (uninitialized) object of class 
                             `LinStatExpectCovar' with p = 1
    *\param expcovinf an initialized object of class `ExpectCovarInfluence'
    *\param splitctrl an object of class `SplitControl'
*/

SEXP R_split(SEXP x, SEXP y, SEXP weights, SEXP orderx, SEXP linexpcov2sample, 
             SEXP expcovinf, SEXP splitctrl) {
             
    SEXP ans, cutpoint, maxstat, statistics;
    
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, cutpoint = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, maxstat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, statistics = allocVector(REALSXP, nrow(x)));
    
    C_split(REAL(x), ncol(x), REAL(y), ncol(y), REAL(weights), nrow(x),
            INTEGER(orderx), splitctrl, linexpcov2sample, expcovinf,
            REAL(cutpoint), REAL(maxstat), REAL(statistics));
    UNPROTECT(1);
    return(ans);
}


/**
    Search for a cutpoint in a unordered factor x maximizing a two-sample
    statistic w.r.t. (the influence function of ) the response variable y.
    *\param codingx the coding of x, i.e. as.numeric(x)
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param codingx the coding of x, i.e. as.numeric(x)
    *\param standstat the vector of the standardized statistics for x, y, 
                      weights 
    *\param splitctrl an object of class `SplitControl'
    *\param linexpcov2sample an (uninitialized) object of class 
                             `LinStatExpectCovar' with p = 1
    *\param expcovinf an initialized object of class `ExpectCovarInfluence'
    *\param cutpoint return value; pointer to a double for the cutpoint in x
    *\param levelset return value; pointer to a p-dim 0/1 integer
    *\param maxstat return value; pointer to a double for the 
                    maximal test statistic
    *\param statistics return value; pointer to a n-dim double for 
                       the statistics
*/

void C_splitcategorical(const int *codingx, int p,
                        const double *y, int q,
                        const double *weights, int n,
                        double *standstat,
                        SEXP splitctrl, SEXP linexpcov2sample, 
                        SEXP expcovinf, double *cutpoint, int *levelset, 
                        double *maxstat, double *statistics) {

    double *tmpx, *tmptmpx, tmp = 0.0;
    int *irank, *ordertmpx, i, j, k, l, jp, chk;

    /* allocate memory */
    tmpx = Calloc(n, double);
    ordertmpx = Calloc(n, int);
    irank = Calloc(p, int);
    tmptmpx = Calloc(n, double);

    /* for all response variables (aka: dummy variables) */
    for (j = 0; j < q; j++) {
    
        jp = j * p;

        /* determine the ranking of the kth level among 
           the standardized statistic: This induced an ordering of the 
           observations */
        for (k = 0; k < p; k++) {
            irank[k] = 1;
            for (l = 0; l < p; l++)
                if (standstat[jp + l] < standstat[jp + k]) irank[k]++;
        }
        
        /* a temporary response variable: the rank of the level */
        for (i = 0; i < n; i++) {
            /* <FIXME> do we have to adjust weights for missing values here??? */
            if (weights[i] == 0.0) {
                tmpx[i] = 0.0;
            } else {
                tmpx[i] = (double) irank[codingx[i] - 1];
            }
            /* </FIXME> */
            tmptmpx[i] = tmpx[i];
            ordertmpx[i] = i + 1;
        }
        
        /* order(dtmpx) */
        rsort_with_index(tmptmpx, ordertmpx, n);

        /* search for a cutpoint (now we do have an ordering) */
        C_split(tmpx, 1, y, q, weights, n, ordertmpx,
                splitctrl, linexpcov2sample,
                expcovinf, cutpoint, maxstat, statistics);

        /* if we have seen an improvement: save this segmentation 
           note: there may be splits with equal goodness */
        chk = 0;
        if (maxstat[0] > tmp) {
            for (k = 0; k < p; k++) {
                if (irank[k] > cutpoint[0]) {
                    levelset[k] = 1;
                    chk += 1;
                } else {
                    levelset[k] = 0;
                }
            }
            tmp = maxstat[0];
        }
        /* <FIXME> make sure that at least one level goes left,
                   C_split may end up with cutpoint > max(irank), why?
           </FIXME>
        */
        /* hm, why did I added 
        if (chk == 0) tmp = 0.0; 
        ??? */
    }
    maxstat[0] = tmp;

    /* free memory */
    Free(tmpx); Free(ordertmpx); Free(irank); Free(tmptmpx);
}


/**
    R-interface to C_splitcategorical (does not handle ordered y's)
    *\param x the values of the x-transformation
    *\param codingx the coding of x, i.e. as.numeric(x)
    *\param y values of the influence function
    *\param weights case weights
    *\param linexpcov2sample an (uninitialized) object of class 
                             `LinStatExpectCovar' with p = 1
    *\param linexpcov an initialized object of class `LinStatExpectCovar'
    *\param expcovinf an initialized object of class `ExpectCovarInfluence'
    *\param splitctrl an object of class `SplitControl'
*/

SEXP R_splitcategorical(SEXP x, SEXP codingx, SEXP y, SEXP weights, 
                        SEXP linexpcov2sample, SEXP linexpcov, 
                        SEXP expcovinf, SEXP splitctrl) {
             
    SEXP ans, cutpoint, maxstat, statistics, levelset;
    double *standstat;

    C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y), REAL(weights), nrow(x),
                    1, GET_SLOT(linexpcov, PL2_expcovinfSym), linexpcov);

    standstat = Calloc(get_dimension(linexpcov), double);
    C_standardize(REAL(GET_SLOT(linexpcov, PL2_linearstatisticSym)),
                  REAL(GET_SLOT(linexpcov, PL2_expectationSym)),
                  REAL(GET_SLOT(linexpcov, PL2_covarianceSym)),
                  get_dimension(linexpcov), get_tol(splitctrl), standstat);

    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, cutpoint = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, maxstat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, statistics = allocVector(REALSXP, nrow(x)));
    SET_VECTOR_ELT(ans, 3, levelset = allocVector(INTSXP, ncol(x)));
    
    C_splitcategorical(INTEGER(codingx), ncol(x), REAL(y), ncol(y), REAL(weights), 
                       nrow(x), standstat, 
                       splitctrl, linexpcov2sample, expcovinf, 
                       REAL(cutpoint), INTEGER(levelset), REAL(maxstat), 
                       REAL(statistics));

    UNPROTECT(1);
    Free(standstat);
    return(ans);
}
