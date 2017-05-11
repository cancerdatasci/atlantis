
/**
    The tree growing recursion
    *\file TreeGrow.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"


int *copy_with_ignoring(int *variables_to_ignore, int ninputs, int index) {
   int i;
   int *new_mask;
   
   if(variables_to_ignore == NULL) {
//     printf("copy_with_ignoring(NULL, %d, %d)\n", ninputs, index);
     return NULL;
   }
   
/*   printf("copy_with_ignoring([");
   for(i = 0;i<ninputs;i++) {
     if(i>0) { 
       printf(", ");
     }
     printf("%d", variables_to_ignore[i]);
   }
   printf("], %d, %d)\n", ninputs, index);
 */
   
   new_mask = Calloc(ninputs, int);
   for(i = 0;i<ninputs;i++) {
     new_mask[i] = i == index ? 1 : variables_to_ignore[i];
   }
   return new_mask;
}

/**
    The main tree growing function, handles the recursion. \n
    *\param node  a list representing the current node
    *\param learnsample an object of class `LearningSample'
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param where a pointer to an integer vector of n-elements
    *\param nodenum a pointer to a integer vector of length 1
    *\param depth an integer giving the depth of the current node
*/

void C_TreeGrow(SEXP node, SEXP learnsample, SEXP fitmem, 
                SEXP controls, int *where, int *nodenum, int depth, int *variables_to_ignore) {

    SEXP weights;
    int nobs, i, stop;
    double *dweights;
    SEXP split;
    int ninputs;
    int *child_variables_to_ignore;
    
    ninputs = get_ninputs(learnsample);
    weights = S3get_nodeweights(node);
    
    /* stop if either stumps have been requested or 
       the maximum depth is exceeded */
    stop = (nodenum[0] == 2 || nodenum[0] == 3) && 
           get_stump(get_tgctrl(controls));
    stop = stop || !check_depth(get_tgctrl(controls), depth);
    
    if (stop)
        C_Node(node, learnsample, weights, fitmem, controls, 1, depth, variables_to_ignore);
    else
        C_Node(node, learnsample, weights, fitmem, controls, 0, depth, variables_to_ignore);
    
    S3set_nodeID(node, nodenum[0]);    
    
    if (!S3get_nodeterminal(node)) {

        C_splitnode(node, learnsample, controls);

        /* determine surrogate splits and split missing values */
        if (get_maxsurrogate(get_splitctrl(controls)) > 0) {
            C_surrogates(node, learnsample, weights, controls, fitmem);
            C_splitsurrogate(node, learnsample);
        }

        split = S3get_primarysplit(node);
        child_variables_to_ignore = copy_with_ignoring(variables_to_ignore, ninputs, S3get_variableID(split)-1);
        nodenum[0] += 1;
        C_TreeGrow(S3get_leftnode(node), learnsample, fitmem, 
                   controls, where, nodenum, depth + 1, child_variables_to_ignore);

        nodenum[0] += 1;                                      
        C_TreeGrow(S3get_rightnode(node), learnsample, fitmem, 
                   controls, where, nodenum, depth + 1, child_variables_to_ignore);
                   
        Free(child_variables_to_ignore);
                   
    } else {
        dweights = REAL(weights);
        nobs = get_nobs(learnsample);
        for (i = 0; i < nobs; i++)
            if (dweights[i] > 0) where[i] = nodenum[0];
    } 
}


/**
    R-interface to C_TreeGrow\n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param where a vector of node indices for each observation
*/

SEXP R_TreeGrow(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls, SEXP where) {
            
     SEXP ans, nweights;
     double *dnweights, *dweights;
     int nobs, i, nodenum = 1;

     GetRNGstate();
     
     nobs = get_nobs(learnsample);
     PROTECT(ans = allocVector(VECSXP, NODE_LENGTH));
     C_init_node(ans, nobs, get_ninputs(learnsample), get_maxsurrogate(get_splitctrl(controls)),
                 ncol(get_predict_trafo(GET_SLOT(learnsample, PL2_responsesSym))));

     nweights = S3get_nodeweights(ans);
     dnweights = REAL(nweights);
     dweights = REAL(weights);
     for (i = 0; i < nobs; i++) dnweights[i] = dweights[i];
     
     C_TreeGrow(ans, learnsample, fitmem, controls, INTEGER(where), &nodenum, 1, NULL);
     
     PutRNGstate();
     
     UNPROTECT(1);
     return(ans);
}
