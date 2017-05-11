
/**
    Some commonly needed utility functions.
    *\file Utils.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"
                
                
/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans) {

    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
                }
            }
        }
    }
}  


/**
    R-interface to C_kronecker
    *\param A matrix
    *\param B matrix
*/
                
SEXP R_kronecker (SEXP A, SEXP B) {

    /*  The Kronecker product, a real (mr x ns) matrix */
    SEXP ans; 
    int *adim, *bdim;

    if (!isReal(A) || !isReal(B)) 
        error("R_kronecker: A and B are not of type REALSXP");

    if (isMatrix(A)) {
        adim = INTEGER(getAttrib(A, R_DimSymbol));
    } else {
        /* assume row vectors */
        adim = Calloc(2, int);
        adim[0] = 1;
        adim[1] = LENGTH(A);
    }
    
    if (isMatrix(B)) {
        bdim = INTEGER(getAttrib(B, R_DimSymbol));
    } else {
        /* assume row vectors */
        bdim = Calloc(2, int);
        bdim[0] = 1;
        bdim[1] = LENGTH(B);
    }

    PROTECT(ans = allocMatrix(REALSXP, 
                              adim[0] * bdim[0], 
                              adim[1] * bdim[1]));
    C_kronecker(REAL(A), adim[0], adim[1], 
                REAL(B), bdim[0], bdim[1], REAL(ans));
    if (!isMatrix(A)) Free(adim); 
    if (!isMatrix(B)) Free(bdim);
    UNPROTECT(1);
    return(ans);
}


/**
    C- and R-interface to La_svd 
    *\param jobu
    *\param jobv
    *\param x
    *\param s
    *\param u
    *\param v
    *\param method
*/

void CR_La_svd(int dim, SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v,
               SEXP method)
{
    int *xdims, n, p, lwork, info = 0;
    int *iwork;
    double *work, *xvals, tmp;
    /* const char * meth; not used*/

    if (!(isString(jobu) && isString(jobv)))
	error(("'jobu' and 'jobv' must be character strings"));
    if (!isString(method))
	error(("'method' must be a character string"));
    /* meth = CHAR(STRING_ELT(method, 0)); not used */
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0]; p = xdims[1];
    xvals = Calloc(n * p, double);
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t) (n * p));

    {
	int ldu = INTEGER(getAttrib(u, R_DimSymbol))[0],
	    ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];
        /* this only works for SQUARE matrices, potentially
           with reduced dimension, so use dim x dim for input and
           output */
        ldu = dim;
        ldvt = dim;
	iwork= (int *) Calloc(8*(n<p ? n : p), int);

	/* ask for optimal size of work array */
	lwork = -1;
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
/* was			 &n, &p, xvals, &n, REAL(s), for the non-square case */
			 &dim, &dim, xvals, &dim, REAL(s),
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 &tmp, &lwork, iwork, &info);
	if (info != 0)
	    error(("error code %d from Lapack routine '%s'"), info, "dgesdd");
	lwork = (int) tmp;
	work = Calloc(lwork, double);
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
			 &dim, &dim, xvals, &dim, REAL(s),
/* was			 &n, &p, xvals, &n, REAL(s), for the non-square case */
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 work, &lwork, iwork, &info);
	if (info != 0)
	    error(("error code %d from Lapack routine '%s'"), info, "dgesdd");
    }
    Free(work); Free(xvals); Free(iwork);
}

/**
    C-interface to CR_La_svd
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

void C_svd (SEXP x, SEXP svdmem) {

    int p, i;
    double *du, *dv;

    if (!isMatrix(x) || !isReal(x))
        error("x is not a real matrix");

    du = REAL(GET_SLOT(svdmem, PL2_uSym));
    dv = REAL(GET_SLOT(svdmem, PL2_vSym));
    p = INTEGER(GET_SLOT(svdmem, PL2_pSym))[0];
    if (nrow(x) < p) error("svd p x error");
    for (i = 0; i < p*p; i++) {
        du[i] = 0.0;
        dv[i] = 0.0;
    }
    CR_La_svd(p, GET_SLOT(svdmem, PL2_jobuSym), 
        GET_SLOT(svdmem, PL2_jobvSym), x, GET_SLOT(svdmem, PL2_sSym), 
        GET_SLOT(svdmem, PL2_uSym), GET_SLOT(svdmem, PL2_vSym), 
        GET_SLOT(svdmem, PL2_methodSym));
    /* return(R_NilValue); */
}


/**
    R-interface to CR_La_svd
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_svd (SEXP x, SEXP svdmem) {

    C_svd(x, svdmem);
    return(R_NilValue);
}


/**
    Reorder the linear statistic, expectation, and covariance
    in a way that elements with zero variance come last. These
    will be ignored in later computation of the MP inverse
    *\param x an object of class `LinStatExpectCovarMPinv'
*/

void C_linexpcovReduce (SEXP x) {

    double *dlinstat, *dexp, *dcov;
    double *dlinstat2, *dexp2, *dcov2;
    int pq, pqn, *zerovar, i, j, itmp, jtmp, sumzv = 0, *dim;

    /* the statistic, expectation, covariance in original dimension */    
    pq = INTEGER(GET_SLOT(x, PL2_dimensionSym))[0];
    dlinstat = REAL(GET_SLOT(x, PL2_linearstatisticSym));
    dexp = REAL(GET_SLOT(x, PL2_expectationSym));
    dcov = REAL(GET_SLOT(x, PL2_covarianceSym));

    /* indicator of zero variance */
    zerovar = Calloc(pq, int);

    /* identify and count zero variances (we can use 0.0 because variances
       corresponding to empty levels are exactly zero */
    for (i = 0; i < pq; i++) {
        if (dcov[i + i * pq] > 0.0) {
            zerovar[i] = 0;
        } else {
            zerovar[i] = 1;
            sumzv++;
        }
    }

    /* if there is any such element and not all variances are zero*/
    if ((sumzv > 0) & (sumzv < pq)) {

        /* do we really need a copy ? */
        pqn = pq - sumzv;
        dlinstat2 = Calloc(pq, double);
        dexp2 = Calloc(pq, double);
        dcov2 = Calloc(pq * pq, double);
        
        /* init */
        for (i = 0; i < pq; i++) {
            dlinstat2[i] = 0.0; 
            dexp2[i] = 0.0; 
            for (j = 0; j < pq; j++)
                dcov2[i + j*pq] = 0.0;
        }
        
        /* overwrite zero variance elements with subsequent elements */
        itmp = 0;
        for (i = 0; i < pq; i++) {
            if (zerovar[i] == 0) {
                dlinstat2[itmp] = dlinstat[i];
                dexp2[itmp] = dexp[i];
                jtmp = 0;
                for (j = 0; j < pq; j++) {
                    if (zerovar[j] == 0) {
                        dcov2[itmp + jtmp * pqn] = dcov[i + j * pq];
                        jtmp++;
                    }
                }
                itmp++;
            }
        }
                                        
        for (i = 0; i < pq; i++) {
            dlinstat[i] = dlinstat2[i]; 
            dexp[i] = dexp2[i]; 
            for (j = 0; j < pq; j++)
                dcov[i + j*pq] = dcov2[i + j * pq];
        }

        /* ATTENTION: This is dangerous but the only way to tell
           svd and friends that we want to use the first pqn elements only! */
        dim = INTEGER(GET_SLOT(x, PL2_dimensionSym));
        dim[0] = pqn;
        /* we reset the original dimension in C_TestStatistic */
        
        Free(dlinstat2);
        Free(dexp2);
        Free(dcov2);
    }
    Free(zerovar);
}


/**
    R-interface to C_linexpcovReduce
    *\param x an object of class `LinStatExpectCovarMPinv'
*/

SEXP R_linexpcovReduce (SEXP x) {

    C_linexpcovReduce(x);
    return(R_NilValue);
}
        

/**
    Moore-Penrose inverse of a matrix
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
    *\param ans return value; an object of class `ExpectCovarMPinv'
*/

void C_MPinv (SEXP x, double tol, SEXP svdmem, SEXP ans) {

    SEXP d, u, vt;
    int i, j, p, pn, k, *positive;
    double *dd, *du, *dvt, *dMPinv;
    double *drank;
    
    drank = REAL(GET_SLOT(ans, PL2_rankSym));
    dMPinv = REAL(GET_SLOT(ans, PL2_MPinvSym));

    C_svd(x, svdmem);
    d = GET_SLOT(svdmem, PL2_sSym);
    dd = REAL(d);
    u = GET_SLOT(svdmem, PL2_uSym);
    du = REAL(u);
    vt = GET_SLOT(svdmem, PL2_vSym);
    dvt = REAL(vt);
    /* this may be the reduced dimension! Use the first p elements only!!!*/
    p = INTEGER(GET_SLOT(svdmem, PL2_pSym))[0];

    if (tol * dd[0] > tol) tol = tol * dd[0];

    positive = Calloc(p, int); 
    
    drank[0] = 0.0;
    for (i = 0; i < p; i++) {
        if (dd[i] > tol) {
            positive[i] = 1;
            drank[0] += 1.0;
        } 
    }
    
    for (j = 0; j < p; j++) {
        if (positive[j]) {
            for (i = 0; i < p; i++)
                du[j * p + i] *= (1 / dd[j]);
        }
    }
    
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++) {
            dMPinv[j * p + i] = 0.0;
            for (k = 0; k < p; k++) {
                if (positive[k])
                    dMPinv[j * p + i] += dvt[i * p + k] * du[p * k + j];
            }
        }
    }

    Free(positive);
}


/**
    R-interface to C_MPinv 
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_MPinv (SEXP x, SEXP tol, SEXP svdmem) {

    SEXP ans;
    int p;

    if (!isMatrix(x) || !isReal(x))
        error("R_MPinv: x is not a real matrix");

    if (nrow(x) != ncol(x)) 
        error("R_MPinv: x is not a square matrix");

    if (!isReal(tol) || LENGTH(tol) != 1)
        error("R_MPinv: tol is not a scalar real");

    p = nrow(x);
    /* potentially, the effective dimension was reduced
    if (p != INTEGER(GET_SLOT(svdmem, PL2_pSym))[0])
        error("R_MPinv: dimensions don't match");
    */

    PROTECT(ans = NEW_OBJECT(MAKE_CLASS("LinStatExpectCovarMPinv")));
    SET_SLOT(ans, PL2_MPinvSym, PROTECT(allocMatrix(REALSXP, p, p)));
    SET_SLOT(ans, PL2_rankSym, PROTECT(allocVector(REALSXP, 1)));
    
    C_MPinv(x, REAL(tol)[0], svdmem, ans);
    
    UNPROTECT(3);
    return(ans);
}


/**
    the maximum of a double vector
    *\param x vector
    *\param n its length
*/


double C_max(const double *x, const int n) {
   double tmp = 0.0;
   int i;
   
   for (i = 0; i < n; i++) {
       if (x[i] > tmp) tmp = x[i];
   }
   return(tmp);
}


/**
    R-interface to C_max
    *\param x numeric vector
*/

SEXP R_max(SEXP x) {

    SEXP ans;
    int n;
    
    if (!isReal(x)) 
        error("R_max: x is not of type REALSXP");
    n = LENGTH(x);
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = C_max(REAL(x), n);
    UNPROTECT(1);
    return(ans);
}


/**
    absolute value 
    *\param x numeric vector
    *\param n length(x)
*/

void C_abs(double *x, int n) {

    int i;
    for (i = 0; i < n; i++) x[i] = fabs(x[i]);
}


/**
    R-interface to C_abs
    *\param x numeric vector
*/

SEXP R_abs(SEXP x) {

    SEXP ans;
    int n;
    
    if (!isReal(x)) 
        error("R_max: x is not of type REALSXP");
    n = LENGTH(x);
    PROTECT(ans = duplicate(x));
    C_abs(REAL(ans), n);
    UNPROTECT(1);
    return(ans);
}


/**
    matrix product x %*% y
    *\param x a matrix
    *\param nrx number of rows of x
    *\param ncx number of cols of x
    *\param y a matrix
    *\param nry number of rows of y
    *\param ncy number of cols of y
    *\param z a matrix of dimension nrx x ncy
*/

void C_matprod(double *x, int nrx, int ncx,
               double *y, int nry, int ncy, double *z)
{
    const char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    int i;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
	                x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}


/**
    R-interface to C_matprod
    *\param x a matrix
    *\param y a matrix
*/

SEXP R_matprod(SEXP x, SEXP y) {

    SEXP ans;
    
    int nrx, ncx, nry, ncy;
    
    nrx = nrow(x);
    ncx = ncol(x);
    nry = nrow(y);
    ncy = ncol(y);

    if (ncx != nry)
        error("R_matprod: dimensions don't match");
    PROTECT(ans = allocMatrix(REALSXP, nrx, ncy));
    C_matprod(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/**
    matrix product x %*% t(y)
    *\param x a matrix
    *\param nrx number of rows of x
    *\param ncx number of cols of x
    *\param y a matrix
    *\param nry number of rows of y
    *\param ncy number of cols of y
    *\param z a matrix of dimension nrx x ncy
*/

void C_matprodT(double *x, int nrx, int ncx,
                double *y, int nry, int ncy, double *z)
{
    const char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    int i;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncy, &one,
	                x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
}


/**
    R-interface to C_matprodT
    *\param x a matrix
    *\param y a matrix
*/

SEXP R_matprodT(SEXP x, SEXP y) {

    SEXP ans;
    int nrx, ncx, nry, ncy;
    
    nrx = nrow(x);
    ncx = ncol(x);
    nry = nrow(y);
    ncy = ncol(y);

    if (ncx != ncy)
        error("R_matprod: dimensions don't match");
    PROTECT(ans = allocMatrix(REALSXP, nrx, nry));
    C_matprodT(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/**
    compute a permutation of a (random subset of) 0:(m-1)
    *\param x an integer vector of length m
    *\param m integer
    *\param k integer
    *\param ans an integer vector of length k
*/    

void C_SampleNoReplace(int *x, int m, int k, int *ans) {
     
    int i, j, n = m;

    for (i = 0; i < m; i++)
        x[i] = i;
    for (i = 0; i < k; i++) {
        j = floor((double) n * unif_rand()); 
        ans[i] = x[j];
        x[j] = x[--n];  
    }
}


/**
    R-interface to C_SampleNoReplace: the permutation case
    *\param m integer
*/    

SEXP R_permute(SEXP m) {
    
    SEXP x, ans;
    int n;
    
    n = INTEGER(m)[0];
    PROTECT(x = allocVector(INTSXP, n));
    PROTECT(ans = allocVector(INTSXP, n));
    C_SampleNoReplace(INTEGER(x), n, n, INTEGER(ans));
    UNPROTECT(2);
    return(ans);
}


/**
    R-interface to C_SampleNoReplace: the subset case
    *\param m integer
    *\param k integer
*/    

SEXP R_rsubset(SEXP m, SEXP k) {
    
    SEXP x, ans;
    int n, j;
    
    n = INTEGER(m)[0];
    j = INTEGER(k)[0];
    PROTECT(x = allocVector(INTSXP, n));
    PROTECT(ans = allocVector(INTSXP, j));
    C_SampleNoReplace(INTEGER(x), n, j, INTEGER(ans));
    UNPROTECT(2);
    return(ans);
}

/* Unequal probability sampling; without-replacement case */

void C_ProbSampleNoReplace(int n, double *p, int *perm,
                           int nans, int *ans)
{
    double rT, mass, totalmass;
    int i, j, k, n1;

    /* Record element identities */
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
	rT = totalmass * unif_rand();
	mass = 0;
	for (j = 0; j < n1; j++) {
	    mass += p[j];
	    if (rT <= mass)
		break;
	}
	ans[i] = perm[j];
	totalmass -= p[j];
	for(k = j; k < n1; k++) {
	    p[k] = p[k + 1];
	    perm[k] = perm[k + 1];
	}
    }
}


/**
    determine if i is element of the integer vector set
    *\param i an integer
    *\param iset a pointer to an integer vector
    *\param p length(iset)
*/

int i_in_set(int i, int *iset, int p) {

    int j, is = 0;
        
    if (p == 0) return(0);
                    
    for (j = 0; j < p; j++) {
        if (iset[j] == i) {  
            is = 1;
            break; 
        }
    }
    return(is);
}

int C_i_in_set(int i, SEXP set) {
    if (LENGTH(set) > 0)
        return(i_in_set(i, INTEGER(set), LENGTH(set)));
    else 
        return(0);
}
    
int nrow(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
}

int ncol(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
}

/* compute index of variable with smallest p-value 
   (and largest test statistic in case two or more p-values coincide -- 
    should not happen anymore since we use 1 - (1 - p)^k for Bonferroni adjustment)
*/
int C_whichmax(double *pvalue, double *teststat, int ninputs) {

    int ans = -1, j;
    double tmppval = 0.0, tmptstat = 0.0;
       
    /* <FIXME> can we switch to the log scale here? </FIXME> */

    tmppval = 0.0;
    tmptstat = 0.0;
    for (j = 0; j < ninputs; j++) {
        if (pvalue[j] > tmppval) {
            ans = j;
            tmppval = pvalue[j];
            tmptstat = teststat[j];
        } else {
            if (pvalue[j] == tmppval && teststat[j] > tmptstat) {  
                ans = j;
                tmppval = pvalue[j];
                tmptstat = teststat[j];
            }
        }
    }
    return(ans);
}

SEXP R_whichmax(SEXP x, SEXP y) {
    SEXP ans;
    
    if (LENGTH(x) != LENGTH(y)) error("different length");
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = C_whichmax(REAL(x), REAL(y), LENGTH(x));
    UNPROTECT(1);
    return(ans);
}

SEXP R_listplus(SEXP a, SEXP b, SEXP which) {

    int na, nb, i, j, *iwhich;
    double *dae, *dbe;
    SEXP ae, be;

    na = LENGTH(a);
    nb = LENGTH(b);
    if (na != nb) error("a and b are of different length");
    
    iwhich = LOGICAL(which);
    
    for (i = 0; i < na; i++) {
        if (iwhich[i]) continue;
        
        ae = VECTOR_ELT(a, i);
        be = VECTOR_ELT(b, i);

        if (LENGTH(ae) != LENGTH(be)) 
            error("elements %d are of different length", i);
            
        if (!isReal(ae) || !isReal(be))
            error("elements %d are not of type double", i);
            
        dae = REAL(ae);
        dbe = REAL(be);
        for (j = 0; j < LENGTH(ae); j++) 
            dae[j] += dbe[j];
    }
    return(a);
}

SEXP R_modify_response(SEXP x, SEXP vf) {

    double *src, *tar;
    int i, n;
    
    src = REAL(x);
    n = LENGTH(x);

    tar = REAL(get_transformation(vf, 1));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_test_trafo(vf));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_predict_trafo(vf));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_variable(vf, 1));
    for (i = 0; i < n; i++)
        tar[i] = src[i];
                                          
    return(R_NilValue);
}

double F77_SUB(unifrnd)(void) { return unif_rand(); }

void C_SampleSplitting(int n, double *prob, int *weights, int k) {

    int i;
    double *tmpprob;
    int *ans, *perm;

    tmpprob = Calloc(n, double);
    perm = Calloc(n, int);
    ans = Calloc(k, int);
    for (i = 0; i < n; i++) tmpprob[i] = prob[i];

    C_ProbSampleNoReplace(n, tmpprob, perm, k, ans);
    for (i = 0; i < n; i++) weights[i] = 0;
    for (i = 0; i < k; i++)
        weights[ans[i] - 1] = 1;
    Free(tmpprob); Free(perm); Free(ans);
}

/**
    Remove weights vector from each node of a tree (in order to save memory)
    \*param subtree a tree
*/ 

void C_remove_weights(SEXP subtree, int removestats) {

    SET_VECTOR_ELT(subtree, S3_WEIGHTS, R_NilValue);
    
    if (!S3get_nodeterminal(subtree)) {
        if (removestats) {
            SET_VECTOR_ELT(VECTOR_ELT(subtree, S3_CRITERION), 
                           S3_iCRITERION, R_NilValue);
            SET_VECTOR_ELT(VECTOR_ELT(subtree, S3_CRITERION), 
                           S3_STATISTICS, R_NilValue);
        }
        C_remove_weights(S3get_leftnode(subtree), removestats);
        C_remove_weights(S3get_rightnode(subtree), removestats);
    }
}

SEXP R_remove_weights(SEXP subtree, SEXP removestats) {

    C_remove_weights(subtree, LOGICAL(removestats)[0]);
    return(R_NilValue);
}

double* C_tempweights(int j, SEXP weights, SEXP fitmem, SEXP inputs) {

    int nobs, *iNAs, i, k;
    double *dw, *dweights;
    SEXP NAs;
    
    dw = REAL(get_weights(fitmem));
    nobs = LENGTH(weights);
    dweights = REAL(weights);
    NAs = get_missings(inputs, j);
    iNAs = INTEGER(NAs);
    if (length(NAs) == 0) return(dw);
    for (i = 0; i < nobs; i++) dw[i] = dweights[i];
    for (k = 0; k < LENGTH(NAs); k++)
        dw[iNAs[k] - 1] = 0.0;
    
    return(dw);
}
