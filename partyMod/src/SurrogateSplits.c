
/**
    Suggorgate splits
    *\file SurrogateSplits.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"

/**
    Search for surrogate splits for bypassing the primary split \n
    *\param node the current node with primary split specified
    *\param learnsample learning sample
    *\param weights the weights associated with the current node
    *\param controls an object of class `TreeControl'
    *\param fitmem an object of class `TreeFitMemory'
    *\todo enable nominal surrogate split variables as well
*/

void C_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, 
                  SEXP fitmem) {

    SEXP x, y, expcovinf; 
    SEXP splitctrl, inputs; 
    SEXP split, thiswhichNA;
    int nobs, ninputs, i, j, k, jselect, maxsurr, *order, nvar = 0;
    double ms, cp, *thisweights, *cutpoint, *maxstat, 
           *splitstat, *dweights, *tweights, *dx, *dy;
    double cut, *twotab, *ytmp, sumw = 0.0;
    
    nobs = get_nobs(learnsample);
    ninputs = get_ninputs(learnsample);
    splitctrl = get_splitctrl(controls);
    maxsurr = get_maxsurrogate(splitctrl);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    jselect = S3get_variableID(S3get_primarysplit(node));
    
    /* (weights > 0) in left node are the new `response' to be approximated */
    y = S3get_nodeweights(VECTOR_ELT(node, S3_LEFT));
    ytmp = Calloc(nobs, double);
    for (i = 0; i < nobs; i++) {
        ytmp[i] = REAL(y)[i];
        if (ytmp[i] > 1.0) ytmp[i] = 1.0;
    }

    for (j = 0; j < ninputs; j++) {
        if (is_nominal(inputs, j + 1)) continue;
        nvar++;
    }
    nvar--;

    if (maxsurr != LENGTH(S3get_surrogatesplits(node)))
        error("nodes does not have %d surrogate splits", maxsurr);
    if (maxsurr > nvar)
        error("cannot set up %d surrogate splits with only %d ordered input variable(s)", 
              maxsurr, nvar);

    tweights = Calloc(nobs, double);
    dweights = REAL(weights);
    for (i = 0; i < nobs; i++) tweights[i] = dweights[i];
    if (has_missings(inputs, jselect)) {
        thiswhichNA = get_missings(inputs, jselect);
        for (k = 0; k < LENGTH(thiswhichNA); k++)
            tweights[INTEGER(thiswhichNA)[k] - 1] = 0.0;
    }

    /* check if sum(weights) > 1 */
    sumw = 0.0;
    for (i = 0; i < nobs; i++) sumw += tweights[i];
    if (sumw < 2.0)
        error("can't implement surrogate splits, not enough observations available");

    expcovinf = GET_SLOT(fitmem, PL2_expcovinfssSym);
    C_ExpectCovarInfluence(ytmp, 1, tweights, nobs, expcovinf);
    
    splitstat = REAL(get_splitstatistics(fitmem));
    /* <FIXME> extend `TreeFitMemory' to those as well ... */
    maxstat = Calloc(ninputs, double);
    cutpoint = Calloc(ninputs, double);
    order = Calloc(ninputs, int);
    /* <FIXME> */
    
    /* this is essentially an exhaustive search */
    /* <FIXME>: we don't want to do this for random forest like trees 
       </FIXME>
     */
    for (j = 0; j < ninputs; j++) {
    
         order[j] = j + 1;
         maxstat[j] = 0.0;
         cutpoint[j] = 0.0;

         /* ordered input variables only (for the moment) */
         if ((j + 1) == jselect || is_nominal(inputs, j + 1))
             continue;

         x = get_variable(inputs, j + 1);

         if (has_missings(inputs, j + 1)) {

             thisweights = C_tempweights(j + 1, weights, fitmem, inputs);

             /* check if sum(weights) > 1 */
             sumw = 0.0;
             for (i = 0; i < nobs; i++) sumw += thisweights[i];
             if (sumw < 2.0) continue;
                 
             C_ExpectCovarInfluence(ytmp, 1, thisweights, nobs, expcovinf);
             
             C_split(REAL(x), 1, ytmp, 1, thisweights, nobs,
                     INTEGER(get_ordering(inputs, j + 1)), splitctrl,
                     GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                     expcovinf, &cp, &ms, splitstat);
         } else {
         
             C_split(REAL(x), 1, ytmp, 1, tweights, nobs,
             INTEGER(get_ordering(inputs, j + 1)), splitctrl,
             GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
             expcovinf, &cp, &ms, splitstat);
         }

         maxstat[j] = -ms;
         cutpoint[j] = cp;
    }


    /* order with respect to maximal statistic */
    rsort_with_index(maxstat, order, ninputs);
    
    twotab = Calloc(4, double);
    
    /* the best `maxsurr' ones are implemented */
    for (j = 0; j < maxsurr; j++) {

        if (is_nominal(inputs, order[j])) continue;
        
        for (i = 0; i < 4; i++) twotab[i] = 0.0;
        cut = cutpoint[order[j] - 1];
        SET_VECTOR_ELT(S3get_surrogatesplits(node), j, 
                       split = allocVector(VECSXP, SPLIT_LENGTH));
        C_init_orderedsplit(split, 0);
        S3set_variableID(split, order[j]);
        REAL(S3get_splitpoint(split))[0] = cut;
        dx = REAL(get_variable(inputs, order[j]));
        dy = REAL(y);

        /* OK, this is a dirty hack: determine if the split 
           goes left or right by the Pearson residual of a 2x2 table.
           I don't want to use the big caliber here 
        */
        for (i = 0; i < nobs; i++) {
            twotab[0] += ((dy[i] == 1) && (dx[i] <= cut)) * tweights[i];
            twotab[1] += (dy[i] == 1) * tweights[i];
            twotab[2] += (dx[i] <= cut) * tweights[i];
            twotab[3] += tweights[i];
        }
        S3set_toleft(split, (int) (twotab[0] - twotab[1] * twotab[2] / 
                     twotab[3]) > 0);
    }
    
    Free(maxstat);
    Free(cutpoint);
    Free(order);
    Free(tweights);
    Free(twotab);
    Free(ytmp);
}

/**
    R-interface to C_surrogates \n
    *\param node the current node with primary split specified
    *\param learnsample learning sample
    *\param weights the weights associated with the current node
    *\param controls an object of class `TreeControl'
    *\param fitmem an object of class `TreeFitMemory'
*/


SEXP R_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, 
                  SEXP fitmem) {

    C_surrogates(node, learnsample, weights, controls, fitmem);
    return(S3get_surrogatesplits(node));
    
}

/**
    Split with missing values \n
    *\param node the current node with primary and surrogate splits 
                 specified
    *\param learnsample learning sample
*/

void C_splitsurrogate(SEXP node, SEXP learnsample) {

    SEXP weights, split, surrsplit;
    SEXP inputs, whichNA, whichNAns;
    double cutpoint, *dx, *dweights, *leftweights, *rightweights;
    int *iwhichNA, k;
    int i, nna, ns;
                    
    weights = S3get_nodeweights(node);
    dweights = REAL(weights);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
            
    leftweights = REAL(S3get_nodeweights(S3get_leftnode(node)));
    rightweights = REAL(S3get_nodeweights(S3get_rightnode(node)));
    surrsplit = S3get_surrogatesplits(node);

    /* if the primary split has any missings */
    split = S3get_primarysplit(node);
    if (has_missings(inputs, S3get_variableID(split))) {

        /* where are the missings? */
        whichNA = get_missings(inputs, S3get_variableID(split));
        iwhichNA = INTEGER(whichNA);
        nna = LENGTH(whichNA);

        /* for all missing values ... */
        for (k = 0; k < nna; k++) {
            ns = 0;
            i = iwhichNA[k] - 1;
            if (dweights[i] == 0) continue;
            
            /* loop over surrogate splits until an appropriate one is found */
            while(TRUE) {
            
                if (ns >= LENGTH(surrsplit)) break;

                split = VECTOR_ELT(surrsplit, ns);
                if (has_missings(inputs, S3get_variableID(split))) {
                    whichNAns = get_missings(inputs, S3get_variableID(split));
                    if (C_i_in_set(i + 1, whichNAns)) {
                        ns++;
                        continue;
                    }
                }

                cutpoint = REAL(S3get_splitpoint(split))[0];
                dx = REAL(get_variable(inputs, S3get_variableID(split)));

                if (S3get_toleft(split)) {
                    if (dx[i] <= cutpoint) {
                        leftweights[i] = dweights[i];
                        rightweights[i] = 0.0;
                    } else {
                        rightweights[i] = dweights[i];
                        leftweights[i] = 0.0;
                    }
                } else {
                    if (dx[i] <= cutpoint) {
                        rightweights[i] = dweights[i];
                        leftweights[i] = 0.0;
                    } else {
                        leftweights[i] = dweights[i];
                        rightweights[i] = 0.0;
                    }
                }
                break;
            }
        }
    }
}
