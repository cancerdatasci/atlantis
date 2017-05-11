
/**
    Node computations
    *\file Node.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"


/**
    Compute prediction of a node
    *\param y the response variable (raw numeric values or dummy encoded factor)
    *\param n number of observations
    *\param q number of columns of y
    *\param weights case weights
    *\param sweights sum of case weights
    *\param ans return value; the q-dimensional predictions
*/
        
void C_prediction(const double *y, int n, int q, const double *weights, 
                  const double sweights, double *ans) {

    int i, j, jn;
    
    for (j = 0; j < q; j++) {
        ans[j] = 0.0;
        jn = j * n;
        for (i = 0; i < n; i++) 
            ans[j] += weights[i] * y[jn + i];
        ans[j] = ans[j] / sweights;
    }
}


void mask_pvalue(double *pvalue, int *variables_to_ignore, int ninputs) {
  int i;
  if(variables_to_ignore == NULL) {
    return;
  }
  
  for(i=0;i<ninputs;i++) {
    if(variables_to_ignore[i]) {
        pvalue[i] = R_NegInf;
    }
  }
}


/**
    The main function for all node computations
    *\param node an initialized node (an S3 object!)
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param TERMINAL logical indicating if this node will
                     be a terminal node
    *\param depth an integer giving the depth of the current node
*/

void C_Node(SEXP node, SEXP learnsample, SEXP weights, 
            SEXP fitmem, SEXP controls, int TERMINAL, int depth, int *variables_to_ignore) {
    
    int nobs, ninputs, jselect, q, j, k, i;
    double mincriterion, sweights, *dprediction;
    double *teststat, *pvalue, smax, cutpoint = 0.0, maxstat = 0.0;
    double *standstat, *splitstat;
    SEXP responses, inputs, x, expcovinf, linexpcov;
    SEXP varctrl, splitctrl, gtctrl, tgctrl, split, testy, predy;
    double *dxtransf, *thisweights;
    int *itable;
    
    nobs = get_nobs(learnsample);
    ninputs = get_ninputs(learnsample);
    varctrl = get_varctrl(controls);
    splitctrl = get_splitctrl(controls);
    gtctrl = get_gtctrl(controls);
    tgctrl = get_tgctrl(controls);
    mincriterion = get_mincriterion(gtctrl);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    testy = get_test_trafo(responses);
    predy = get_predict_trafo(responses);
    q = ncol(testy);

    /* <FIXME> we compute C_GlobalTest even for TERMINAL nodes! </FIXME> */

    /* compute the test statistics and the node criteria for each input */        
    C_GlobalTest(learnsample, weights, fitmem, varctrl,
                 gtctrl, get_minsplit(splitctrl), 
                 REAL(S3get_teststat(node)), REAL(S3get_criterion(node)), depth);
    
    /* sum of weights: C_GlobalTest did nothing if sweights < mincriterion */
    sweights = REAL(GET_SLOT(GET_SLOT(fitmem, PL2_expcovinfSym), 
                             PL2_sumweightsSym))[0];
    REAL(VECTOR_ELT(node, S3_SUMWEIGHTS))[0] = sweights;

    /* compute the prediction of this node */
    dprediction = REAL(S3get_prediction(node));

    /* <FIXME> feed raw numeric values OR dummy encoded factors as y 
       Problem: what happens for survival times ? */
    C_prediction(REAL(predy), nobs, ncol(predy), REAL(weights), 
                     sweights, dprediction);
    /* </FIXME> */

    teststat = REAL(S3get_teststat(node));
    pvalue = REAL(S3get_criterion(node));

    mask_pvalue(pvalue, variables_to_ignore, ninputs);

    /* try the two out of ninputs best inputs variables */
    /* <FIXME> be more flexible and add a parameter controlling
               the number of inputs tried </FIXME> */
    for (j = 0; j < 2; j++) {

        smax = C_max(pvalue, ninputs);
        REAL(S3get_maxcriterion(node))[0] = smax;
    
        /* if the global null hypothesis was rejected */
        if (smax > mincriterion && !TERMINAL) {

            /* the input variable with largest association to the response */
            jselect = C_whichmax(pvalue, teststat, ninputs) + 1;

            /* get the raw numeric values or the codings of a factor */
            x = get_variable(inputs, jselect);
            if (has_missings(inputs, jselect)) {
                expcovinf = GET_SLOT(get_varmemory(fitmem, jselect), 
                                    PL2_expcovinfSym);
                thisweights = C_tempweights(jselect, weights, fitmem, inputs);
            } else {
                expcovinf = GET_SLOT(fitmem, PL2_expcovinfSym);
                thisweights = REAL(weights);
            }

            /* <FIXME> handle ordered factors separatly??? </FIXME> */
            if (!is_nominal(inputs, jselect)) {
            
                /* search for a split in a ordered variable x */
                split = S3get_primarysplit(node);
                
                /* check if the n-vector of splitstatistics 
                   should be returned for each primary split */
                if (get_savesplitstats(tgctrl)) {
                    C_init_orderedsplit(split, nobs);
                    splitstat = REAL(S3get_splitstatistics(split));
                } else {
                    C_init_orderedsplit(split, 0);
                    splitstat = REAL(get_splitstatistics(fitmem));
                }

                C_split(REAL(x), 1, REAL(testy), q, thisweights, nobs,
                        INTEGER(get_ordering(inputs, jselect)), splitctrl, 
                        GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                        expcovinf, REAL(S3get_splitpoint(split)), &maxstat,
                        splitstat);
                S3set_variableID(split, jselect);
             } else {
           
                 /* search of a set of levels (split) in a numeric variable x */
                 split = S3get_primarysplit(node);
                 
                /* check if the n-vector of splitstatistics 
                   should be returned for each primary split */
                if (get_savesplitstats(tgctrl)) {
                    C_init_nominalsplit(split, 
                        LENGTH(get_levels(inputs, jselect)), 
                        nobs);
                    splitstat = REAL(S3get_splitstatistics(split));
                } else {
                    C_init_nominalsplit(split, 
                        LENGTH(get_levels(inputs, jselect)), 
                        0);
                    splitstat = REAL(get_splitstatistics(fitmem));
                }
          
                 linexpcov = get_varmemory(fitmem, jselect);
                 standstat = Calloc(get_dimension(linexpcov), double);
                 C_standardize(REAL(GET_SLOT(linexpcov, 
                                             PL2_linearstatisticSym)),
                               REAL(GET_SLOT(linexpcov, PL2_expectationSym)),
                               REAL(GET_SLOT(linexpcov, PL2_covarianceSym)),
                               get_dimension(linexpcov), get_tol(splitctrl), 
                               standstat);
 
                 C_splitcategorical(INTEGER(x), 
                                    LENGTH(get_levels(inputs, jselect)), 
                                    REAL(testy), q, thisweights, 
                                    nobs, standstat, splitctrl, 
                                    GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                                    expcovinf, &cutpoint, 
                                    INTEGER(S3get_splitpoint(split)),
                                    &maxstat, splitstat);

                 /* compute which levels of a factor are available in this node 
                    (for printing) later on. A real `table' for this node would
                    induce too much overhead here. Maybe later. */
                    
                 itable = INTEGER(S3get_table(split));
                 dxtransf = REAL(get_transformation(inputs, jselect));
                 for (k = 0; k < LENGTH(get_levels(inputs, jselect)); k++) {
                     itable[k] = 0;
                     for (i = 0; i < nobs; i++) {
                         if (dxtransf[k * nobs + i] * thisweights[i] > 0) {
                             itable[k] = 1;
                             continue;
                         }
                     }
                 }

                 Free(standstat);
            }
            if (maxstat == 0) {
                if (j == 1) {          
                    S3set_nodeterminal(node);
                } else {
                    /* do not look at jselect in next iteration */
                    pvalue[jselect - 1] = R_NegInf;
                }
            } else {
                S3set_variableID(split, jselect);
                break;
            }
        } else {
            S3set_nodeterminal(node);
            break;
        }
    }
}       


/**
    R-interface to C_Node
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
*/

SEXP R_Node(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls) {
            
     SEXP ans;
     
     PROTECT(ans = allocVector(VECSXP, NODE_LENGTH));
     C_init_node(ans, get_nobs(learnsample), get_ninputs(learnsample), 
                 get_maxsurrogate(get_splitctrl(controls)),
                 ncol(get_predict_trafo(GET_SLOT(learnsample, PL2_responsesSym))));

     C_Node(ans, learnsample, weights, fitmem, controls, 0, 1, NULL);
     UNPROTECT(1);
     return(ans);
}
