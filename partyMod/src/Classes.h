
extern SEXP
    PL2_expectationSym,
    PL2_covarianceSym,
    PL2_linearstatisticSym,
    PL2_expcovinfSym,
    PL2_expcovinfssSym,
    PL2_sumweightsSym,
    PL2_dimensionSym,
    PL2_MPinvSym,
    PL2_rankSym,
    PL2_svdmemSym,
    PL2_methodSym,
    PL2_jobuSym, 
    PL2_jobvSym, 
    PL2_uSym,
    PL2_vSym,
    PL2_sSym,
    PL2_pSym,
    PL2_teststatSym,
    PL2_pvalueSym,
    PL2_tolSym,
    PL2_maxptsSym,
    PL2_absepsSym,
    PL2_relepsSym,
    PL2_minsplitSym,
    PL2_minbucketSym,
    PL2_minprobSym,
    PL2_variablesSym, 
    PL2_transformationsSym, 
    PL2_is_nominalSym, 
    PL2_is_ordinalSym, 
    PL2_is_censoredSym, 
    PL2_orderingSym, 
    PL2_levelsSym, 
    PL2_scoresSym, 
    PL2_has_missingsSym, 
    PL2_whichNASym, 
    PL2_test_trafoSym, 
    PL2_predict_trafoSym, 
    PL2_nobsSym, 
    PL2_ninputsSym,
    PL2_linexpcov2sampleSym, 
    PL2_weightsSym, 
    PL2_varmemorySym,
    PL2_linexpcov2sampleSym, 
    PL2_weightsSym, 
    PL2_varmemorySym, 
    PL2_responsesSym, 
    PL2_inputsSym,
    PL2_testtypeSym, 
    PL2_nresampleSym,
    PL2_varctrlSym, 
    PL2_splitctrlSym, 
    PL2_gtctrlSym,
    PL2_mincriterionSym,
    PL2_randomsplitsSym,
    PL2_mtrySym,
    PL2_dontuseSym,
    PL2_dontusetmpSym,
    PL2_stumpSym,
    PL2_tgctrlSym,
    PL2_ntreeSym,
    PL2_replaceSym,
    PL2_fractionSym,
    PL2_traceSym,
    PL2_dropcriterionSym,
    PL2_compressSym,
    PL2_expandSym,
    PL2_varOnceSym;
            
int get_dimension(SEXP object);
int get_teststat(SEXP object);
double get_tol(SEXP object);
int get_pvalue(SEXP object);
int get_maxpts(SEXP object);
double get_abseps(SEXP object);
double get_releps(SEXP object);
double get_minsplit(SEXP object);
double get_minprob(SEXP object);
double get_minbucket(SEXP object);
SEXP get_transformation(SEXP object, int variable);
SEXP get_test_trafo(SEXP object);
SEXP get_predict_trafo(SEXP object);
SEXP get_variable(SEXP object, int variable);
int is_nominal(SEXP object, int variable);
int is_ordinal(SEXP object, int variable);
int is_censored(SEXP object, int variable);
int has_missings(SEXP object, int variable);
SEXP get_missings(SEXP object, int variable);
SEXP get_ordering(SEXP object, int variable);
SEXP get_levels(SEXP object, int variable);
SEXP get_scores(SEXP object, int variable); 
SEXP get_whichNA(SEXP object, int variable);
SEXP get_varmemory(SEXP object, int variable);
int get_nobs(SEXP object);
int get_ninputs(SEXP object);
SEXP get_weights(SEXP object);
int get_testtype(SEXP object);
int get_nresample(SEXP object);
SEXP get_varctrl(SEXP object);
SEXP get_splitctrl(SEXP object);
SEXP get_gtctrl(SEXP object);
double get_mincriterion(SEXP object);
int get_randomsplits(SEXP object);
SEXP get_mtry(SEXP object);
SEXP get_dontuse(SEXP object);
SEXP get_dontusetmp(SEXP object);
int get_stump(SEXP object);
int get_maxsurrogate(SEXP object);
SEXP get_tgctrl(SEXP object);
SEXP get_splitstatistics(SEXP object);
int get_savesplitstats(SEXP object);
int check_depth(SEXP object, int depth);
int get_ntree(SEXP object);
int get_replace(SEXP object);
double get_fraction(SEXP object);
int get_trace(SEXP object);
int get_dropcriterion(SEXP object);
SEXP get_compress(SEXP object);
SEXP get_expand(SEXP object);
int get_only_use_variable_once(SEXP object);
