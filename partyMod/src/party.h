
/* the global header file for the `party' package */

/* include R header files */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */

/* include private header files: this need to be restricted */

#include "Classes.h"
#include "Utils.h"
#include "mvt.h"
#include "LinearStatistic.h"
#include "TestStatistic.h"
#include "Distributions.h"
#include "Convenience.h"
#include "S3Classes.h"
#include "IndependenceTest.h"
#include "Splits.h"
#include "Node.h"
#include "Predict.h"
#include "SurrogateSplits.h"
#include "TreeGrow.h"

/* constants, basically the length of lists representing S3 classes
   and the position of certain elements  */

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

/* type of test statistic */
#define MAXABS			1
#define QUADFORM		2

/* type of criterion to be _maximized_! */
#define BONFERRONI		1
#define MONTECARLO		2
#define AGGREGATED		3
#define UNIVARIATE		4
#define TESTSTATISTIC		5
