
void C_GlobalTest(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP varctrl, 
                  SEXP gtestctrl, double minsplit, double *teststat, double *criterion, int depth);
void C_TeststatPvalue(const SEXP linexpcov, const SEXP varctrl,
                      double *ans_teststat, double *ans_pvalue);
void C_TeststatCriterion(const SEXP linexpcov, const SEXP varctrl, 
                         double *ans_teststat, double *ans_criterion);
