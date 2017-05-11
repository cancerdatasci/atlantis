
/**
    S3 classes for dealing with nodes and splits
    *\file S3Classes.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"
                
void C_init_node(SEXP node, int nobs, int ninputs, int nsurr, int q) {

    SEXP nodeID, weights, criterion, primarysplit, surrogatesplits, 
         terminal, prediction;

    if (LENGTH(node) < NODE_LENGTH)
        error("node is not a list with at least %d elements", NODE_LENGTH);
        
    SET_VECTOR_ELT(node, S3_NODEID, nodeID = allocVector(INTSXP, 1));
    if (nobs > 0) 
        SET_VECTOR_ELT(node, S3_WEIGHTS, weights = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(node, S3_WEIGHTS, R_NilValue);
    SET_VECTOR_ELT(node, S3_SUMWEIGHTS, allocVector(REALSXP, 1));
    SET_VECTOR_ELT(node, S3_CRITERION, 
        criterion = allocVector(VECSXP, CRITERION_LENGTH));
    /* teststats */
    SET_VECTOR_ELT(criterion, S3_STATISTICS, allocVector(REALSXP, ninputs)); 
    /* criterion, aka pvalues */
    SET_VECTOR_ELT(criterion, S3_iCRITERION, allocVector(REALSXP, ninputs));
    /* max(criterion) */
    SET_VECTOR_ELT(criterion, S3_MAXCRITERION, allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(node, S3_TERMINAL, terminal = allocVector(LGLSXP, 1));
    INTEGER(terminal)[0] = 0;
    SET_VECTOR_ELT(node, S3_PSPLIT, 
        primarysplit = allocVector(VECSXP, SPLIT_LENGTH));
    SET_VECTOR_ELT(node, S3_SSPLIT, 
                   surrogatesplits = allocVector(VECSXP, nsurr));
    SET_VECTOR_ELT(node, S3_PREDICTION, prediction = allocVector(REALSXP, q));

}

void S3set_nodeID(SEXP node, int nodeID) {
    INTEGER(VECTOR_ELT(node, S3_NODEID))[0] = nodeID;
}

int S3get_nodeID(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, S3_NODEID))[0]);
}

SEXP S3get_nodeweights(SEXP node) {
    SEXP ans;
    
    ans = VECTOR_ELT(node, S3_WEIGHTS);
    if (ans == R_NilValue)
        error("node has no weights element"); 
    return(VECTOR_ELT(node, S3_WEIGHTS));
}

double S3get_sumweights(SEXP node) {
    return(REAL(VECTOR_ELT(node, S3_SUMWEIGHTS))[0]);
}

SEXP S3get_teststat(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_STATISTICS));
}

SEXP S3get_criterion(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_iCRITERION));
}

SEXP S3get_maxcriterion(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_MAXCRITERION));
}

void S3set_nodeterminal(SEXP node) {
    INTEGER(VECTOR_ELT(node, S3_TERMINAL))[0] = 1;
}

int S3get_nodeterminal(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, S3_TERMINAL))[0]);
}

SEXP S3get_primarysplit(SEXP node) {
    return(VECTOR_ELT(node, S3_PSPLIT));
}

SEXP S3get_surrogatesplits(SEXP node) {
    return(VECTOR_ELT(node, S3_SSPLIT));
}

SEXP S3get_prediction(SEXP node) {
    return(VECTOR_ELT(node, S3_PREDICTION));
}

SEXP S3get_leftnode(SEXP node) {
    return(VECTOR_ELT(node, S3_LEFT));
}

SEXP S3get_rightnode(SEXP node) {
    return(VECTOR_ELT(node, S3_RIGHT));
}

void C_init_orderedsplit(SEXP split, int nobs) {
    
    SEXP variableID, splitpoint, splitstatistics, ordered, toleft;
    
    if (LENGTH(split) < SPLIT_LENGTH)
        error("split is not a list with at least %d elements", SPLIT_LENGTH);
        
    SET_VECTOR_ELT(split, S3_VARIABLEID, 
                   variableID = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(split, S3_ORDERED, 
                    ordered = allocVector(LGLSXP, 1));
    INTEGER(ordered)[0] = 1;
    SET_VECTOR_ELT(split, S3_SPLITPOINT, 
                   splitpoint = allocVector(REALSXP, 1));
    if (nobs > 0)
        SET_VECTOR_ELT(split, S3_SPLITSTATISTICS, 
                       splitstatistics = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(split, S3_SPLITSTATISTICS, R_NilValue);
    SET_VECTOR_ELT(split, S3_TOLEFT, toleft = allocVector(INTSXP, 1));
    INTEGER(toleft)[0] = 1;
    SET_VECTOR_ELT(split, S3_TABLE, R_NilValue);
}

void C_init_nominalsplit(SEXP split, int nlevels, int nobs) {
    
    SEXP variableID, splitpoint, splitstatistics, ordered, toleft, table;
    
    if (LENGTH(split) < SPLIT_LENGTH)
        error("split is not a list with at least %d elements", SPLIT_LENGTH);

    SET_VECTOR_ELT(split, S3_VARIABLEID, variableID = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(split, S3_ORDERED, ordered = allocVector(LGLSXP, 1));
    INTEGER(ordered)[0] = 0;
    SET_VECTOR_ELT(split, S3_SPLITPOINT, 
                   splitpoint = allocVector(INTSXP, nlevels));
    if (nobs > 0)
        SET_VECTOR_ELT(split, S3_SPLITSTATISTICS, 
                       splitstatistics = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(split, S3_SPLITSTATISTICS, R_NilValue);
    SET_VECTOR_ELT(split, S3_TOLEFT, toleft = allocVector(INTSXP, 1));
    INTEGER(toleft)[0] = 1;
    SET_VECTOR_ELT(split, S3_TABLE, table = allocVector(INTSXP, nlevels));
}

void S3set_variableID(SEXP split, int variableID) {
    INTEGER(VECTOR_ELT(split, S3_VARIABLEID))[0] = variableID;
}

int S3get_variableID(SEXP split) {
    return(INTEGER(VECTOR_ELT(split, S3_VARIABLEID))[0]);
}

int S3is_ordered(SEXP split) {
    return(INTEGER(VECTOR_ELT(split, S3_ORDERED))[0]);
}

void S3set_ordered(SEXP split) {
    INTEGER(VECTOR_ELT(split, S3_ORDERED))[0] = 1;
}

void S3set_nominal(SEXP split) {
    INTEGER(VECTOR_ELT(split, S3_ORDERED))[0] = 0;
}

int S3get_toleft(SEXP split) {
    return(INTEGER(VECTOR_ELT(split, S3_TOLEFT))[0]);
}

void S3set_toleft(SEXP split, int left) {
    /* <FIXME> use LOGICAL here? </FIXME> */
    INTEGER(VECTOR_ELT(split, S3_TOLEFT))[0] = left;
}

SEXP S3get_splitpoint(SEXP split) {
   return(VECTOR_ELT(split, S3_SPLITPOINT));
}
   
SEXP S3get_splitstatistics(SEXP split) {
   SEXP ans;
   
   ans = VECTOR_ELT(split, S3_SPLITSTATISTICS);
   if (ans == R_NilValue)
       error("split does not have a splitstatistics element");
   return(ans);
}

SEXP S3get_table(SEXP split) {
   SEXP ans;
   
   ans = VECTOR_ELT(split, S3_TABLE);
   if (ans == R_NilValue)
       error("split does not have a table element");
   return(ans);
}
