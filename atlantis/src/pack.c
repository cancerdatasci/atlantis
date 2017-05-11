#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdlib.h>

/* S3 list elements in `splittingNode's */
#define S3_NODEID		0    /* nodeID */
#define S3_WEIGHTS		1    /* weights */
#define S3_CRITERION		2    /* criterion */
#define S3_TERMINAL		3    /* terminal */
#define S3_PSPLIT		4    /* psplit */
#define S3_SSPLIT		5    /* ssplit */
#define S3_PREDICTION		6    /* prediction */
#define S3_LEFT			7    /* left */
#define S3_RIGHT		8    /* right */
#define S3_SUMWEIGHTS           9    /* sum of weights in this node */
#define NODE_LENGTH		10   /* 9 elements in total */

/* S3 list elements in `criterion' element of `SplittingNode's */
#define S3_STATISTICS		0    /* statistics */
#define S3_iCRITERION		1    /* criterion */
#define S3_MAXCRITERION		2    /* max(criterion) */
#define CRITERION_LENGTH	3    /* 3 elements in total */

/* S3 list elements in `orderedSplit's or `nominalSplit's */
#define S3_VARIABLEID		0    /* variableID */
#define S3_ORDERED		1    /* ordered */
#define S3_SPLITPOINT		2    /* splitpoint */
#define S3_SPLITSTATISTICS	3    /* splitstatistics */
#define S3_TOLEFT		4    /* toleft */
#define S3_TABLE                5    /* table for nominal splits */
#define SPLIT_LENGTH		6    /* 6 elements in total */

typedef struct PackedReals {
  int length;
  double values[1];
} PackedReals;

typedef struct Split {
  int variableid;
  int ordered;
  double splitpoint;
  int toleft;
} Split;

typedef struct PackedSplittingNode {
  int nodeid;
  int terminal;
  Split psplit;
  Split ssplit;
  PackedReals *statistics;
  PackedReals *icriterion;
  double maxcriterion;
  double prediction;
  struct PackedSplittingNode *left;
  struct PackedSplittingNode *right;
  double sumweights;
} PackedSplittingNode;


static long bytes_allocated = 0;
static long nodes_allocated = 0;
static long nodes_freed = 0;
static long pack_calls = 0;

void packSplit(Split *s, SEXP ss) {
  //printf("packSplit\n");
  //printf("ss = %p (%d) ==? %p\n", ss, 0, R_NilValue);
  if(ss == R_NilValue || VECTOR_ELT(ss, S3_VARIABLEID) == R_NilValue) {
    s->variableid = -1;
    return;
  }
  
  s->variableid = INTEGER(VECTOR_ELT(ss, S3_VARIABLEID))[0];
//  printf("packSplit 1\n");
  s->ordered = INTEGER(VECTOR_ELT(ss, S3_ORDERED))[0];
//  printf("packSplit 2\n");
  s->splitpoint = REAL(VECTOR_ELT(ss, S3_SPLITPOINT))[0];
//  printf("packSplit 3\n");
  s->toleft = INTEGER(VECTOR_ELT(ss, S3_TOLEFT))[0];
//  printf("packSplit 4\n");
//  S3_SPLITSTATISTICS R_NilValue
//  S3_TABLE
}

PackedReals *pack_reals(SEXP real_vect) {
  int i;
  int len;
  PackedReals *packed;
  double *dst;
  double *src;
  
  if(real_vect == R_NilValue) { 
    return NULL;
  }
  
  len =  length(real_vect);
  packed  = calloc(sizeof(PackedReals) + sizeof(double) * len, 1);
  packed->length = len;
  dst = packed->values;
  src = REAL(real_vect);
  for(i=0;i<len;i++) {
    *dst = *src;
    dst++;
    src++;
  }
  
  return packed;
}

SEXP unpack_reals(PackedReals *packed) {
  SEXP ans;
  double *dst;
  double *src;
  int i, len;
  
  if(packed == NULL) {
    return R_NilValue;
  }
  
  src = packed->values;

  PROTECT(ans = allocVector(REALSXP, packed->length));
  dst = REAL(ans);
  len = packed->length;
  for(i=0;i<len;i++) {
    *dst = *src;
    dst++;
    src++;
  }
  
  UNPROTECT(1);

  return ans;
}

SEXP R_print_pack_alloc_summary(void) {
  printf("bytes_allocated=%ld, nodes_allocated=%ld, nodes_freed=%ld, pack_calls=%ld\n", bytes_allocated, nodes_allocated, nodes_freed, pack_calls);
  return R_NilValue;
}

PackedSplittingNode *packNode(SEXP node) {
  pack_calls ++;
  
  //printf("packNode = %p, %p\n", node, R_NilValue);
  if(node == R_NilValue) {
    return NULL;
  }

  nodes_allocated ++;
  
  PackedSplittingNode * p = calloc(1,sizeof(PackedSplittingNode));
  bytes_allocated += sizeof(PackedSplittingNode);

  //printf("1\n");
  p->nodeid = INTEGER(VECTOR_ELT(node, S3_NODEID))[0];
  //printf("2\n");
  p->terminal = INTEGER(VECTOR_ELT(node, S3_TERMINAL))[0];
  //printf("3\n");
  packSplit(&(p->psplit), VECTOR_ELT(node, S3_PSPLIT));
//  packSplit(&(p->ssplit), VECTOR_ELT(node, S3_SSPLIT));
//  printf("4\n");

// Are these needed?
//  p->statistics = pack_reals(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_STATISTICS));
//  p->icriterion = pack_reals(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_iCRITERION));

  p->statistics = NULL;
  p->icriterion = NULL;
  p->maxcriterion = REAL(VECTOR_ELT(VECTOR_ELT(node, S3_CRITERION), S3_MAXCRITERION))[0];
  p->prediction = REAL(VECTOR_ELT(node, S3_PREDICTION))[0];
  //printf("index %d\n", S3_LEFT);
  p->left = packNode(VECTOR_ELT(node, S3_LEFT));
//  printf("5\n");
  //printf("index %d\n", S3_RIGHT);
  p->right = packNode(VECTOR_ELT(node, S3_RIGHT));
//  printf("6\n");
  p->sumweights = REAL(VECTOR_ELT(node, S3_SUMWEIGHTS))[0];
//  printf("7\n");
  //printf("up\n");  
  
  return p;
}

void freeNode(PackedSplittingNode *p) {
//  printf("Freeing %p\n", p);
  if(p->left != NULL) {
    freeNode(p->left);
  }
  if(p->right != NULL) {
    freeNode(p->right);
  }
  if(p->statistics != NULL) {
    free(p->statistics);
  }
  if(p->icriterion != NULL) {
    free(p->icriterion);
  }
  bytes_allocated -= sizeof(PackedSplittingNode);
  nodes_freed ++;
  free(p);
}

void freePackedSplittingNode(SEXP ptr) {
  if(!R_ExternalPtrAddr(ptr)) return;

  PackedSplittingNode *p = (PackedSplittingNode *) R_ExternalPtrAddr(ptr);
  freeNode(p);
  R_ClearExternalPtr(ptr);
}

SEXP R_packNode(SEXP exp) {
  SEXP ptr;
  
  PackedSplittingNode *packedNode = packNode(exp);
  ptr = R_MakeExternalPtr(packedNode, install("packedSplittingNode"), R_NilValue);
  PROTECT(ptr);
  R_RegisterCFinalizerEx(ptr, freePackedSplittingNode, FALSE);
  UNPROTECT(1);

  return ptr;
}

SEXP unpackSplit(Split *split) {
  SEXP ans, variableid, ordered, splitpoint, toleft;
  
  if(split == NULL || split->variableid == -1) {
//    return R_NilValue;
    return allocVector(VECSXP, SPLIT_LENGTH);
  }
  
  PROTECT(ans = allocVector(VECSXP, SPLIT_LENGTH));

  SET_VECTOR_ELT(ans, S3_VARIABLEID, variableid = allocVector(INTSXP,1));
  INTEGER(variableid)[0] = split->variableid;
  
  SET_VECTOR_ELT(ans, S3_ORDERED, ordered = allocVector(LGLSXP, 1));
  INTEGER(ordered)[0] = split->ordered;
  
  SET_VECTOR_ELT(ans, S3_SPLITPOINT, splitpoint = allocVector(REALSXP, 1));
  REAL(splitpoint)[0] = split->splitpoint;
  
  SET_VECTOR_ELT(ans, S3_SPLITSTATISTICS, R_NilValue);
  
  SET_VECTOR_ELT(ans, S3_TOLEFT, toleft = allocVector(INTSXP, 1));
  INTEGER(toleft)[0] = split->toleft;
  
  SET_VECTOR_ELT(ans, S3_TABLE, R_NilValue);

  UNPROTECT(1);
  return ans;
}


SEXP unpackNode(PackedSplittingNode *node, int first) {
  SEXP ans, nodeid, terminal, prediction, sumweights, criterion, maxcriterion;
  if(node == NULL) return R_NilValue;
  
  if(first) {
    PROTECT(ans = allocVector(VECSXP, NODE_LENGTH+1));
  } else {
    PROTECT(ans = allocVector(VECSXP, NODE_LENGTH));
  }
  //printf("unpack1\n");  
  SET_VECTOR_ELT(ans, S3_NODEID, nodeid = allocVector(INTSXP,1));
  INTEGER(nodeid)[0] = node->nodeid;
  
  //printf("unpack2\n");  
  SET_VECTOR_ELT(ans, S3_TERMINAL, terminal = allocVector(LGLSXP, 1));
  INTEGER(terminal)[0] = node->terminal;

  SET_VECTOR_ELT(ans, S3_CRITERION, criterion = allocVector(VECSXP, CRITERION_LENGTH));
  SET_VECTOR_ELT(criterion, S3_MAXCRITERION, maxcriterion = allocVector(REALSXP, 1));
  REAL(maxcriterion)[0] = node->maxcriterion;
  SET_VECTOR_ELT(criterion, S3_STATISTICS, unpack_reals(node->statistics));
  SET_VECTOR_ELT(criterion, S3_iCRITERION, unpack_reals(node->icriterion));
  
  //printf("unpack3\n");  
  SET_VECTOR_ELT(ans, S3_PSPLIT, unpackSplit(&(node->psplit)));
  //SET_VECTOR_ELT(ans, S3_SSPLIT, unpackSplit(&(node->ssplit)));
  SET_VECTOR_ELT(ans, S3_SSPLIT, allocVector(VECSXP, 0));
  
  //printf("unpack4\n");  
  SET_VECTOR_ELT(ans, S3_PREDICTION, prediction = allocVector(REALSXP, 1));
  REAL(prediction)[0] = node->prediction;
  
  //printf("unpack5\n");  
  SET_VECTOR_ELT(ans, S3_LEFT, unpackNode(node->left, 0));
  //printf("unpack6\n");  
  SET_VECTOR_ELT(ans, S3_RIGHT, unpackNode(node->right, 0));

  SET_VECTOR_ELT(ans, S3_SUMWEIGHTS, sumweights = allocVector(REALSXP, 1));
  REAL(sumweights)[0] = node->sumweights;
  
  //printf("unpack2\n");  
  UNPROTECT(1);
  return ans;
}

SEXP R_unpackNode(SEXP exp) {
  PackedSplittingNode *p = (PackedSplittingNode *) R_ExternalPtrAddr(exp);
  return unpackNode(p, 1);
}
