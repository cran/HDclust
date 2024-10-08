// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_findModes
S4 rcpp_findModes(NumericMatrix dataTranspose, S4 HmmVb, IntegerVector nthread);
RcppExport SEXP _HDclust_rcpp_findModes(SEXP dataTransposeSEXP, SEXP HmmVbSEXP, SEXP nthreadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dataTranspose(dataTransposeSEXP);
    Rcpp::traits::input_parameter< S4 >::type HmmVb(HmmVbSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nthread(nthreadSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_findModes(dataTranspose, HmmVb, nthread));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_clust
S4 rcpp_clust(NumericMatrix dataTranspose, S4 HmmVb, Nullable<List> rfsClust_, List control, IntegerVector nthread);
RcppExport SEXP _HDclust_rcpp_clust(SEXP dataTransposeSEXP, SEXP HmmVbSEXP, SEXP rfsClust_SEXP, SEXP controlSEXP, SEXP nthreadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dataTranspose(dataTransposeSEXP);
    Rcpp::traits::input_parameter< S4 >::type HmmVb(HmmVbSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type rfsClust_(rfsClust_SEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nthread(nthreadSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_clust(dataTranspose, HmmVb, rfsClust_, control, nthread));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_trainHmmVb
S4 rcpp_trainHmmVb(NumericMatrix dataTranspose, const RObject& VbStructure, const List& searchControl, const List& trainControl, IntegerVector nthread, Function VB, Function HMM, Function HMMVB, bool bprint);
RcppExport SEXP _HDclust_rcpp_trainHmmVb(SEXP dataTransposeSEXP, SEXP VbStructureSEXP, SEXP searchControlSEXP, SEXP trainControlSEXP, SEXP nthreadSEXP, SEXP VBSEXP, SEXP HMMSEXP, SEXP HMMVBSEXP, SEXP bprintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dataTranspose(dataTransposeSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type VbStructure(VbStructureSEXP);
    Rcpp::traits::input_parameter< const List& >::type searchControl(searchControlSEXP);
    Rcpp::traits::input_parameter< const List& >::type trainControl(trainControlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nthread(nthreadSEXP);
    Rcpp::traits::input_parameter< Function >::type VB(VBSEXP);
    Rcpp::traits::input_parameter< Function >::type HMM(HMMSEXP);
    Rcpp::traits::input_parameter< Function >::type HMMVB(HMMVBSEXP);
    Rcpp::traits::input_parameter< bool >::type bprint(bprintSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_trainHmmVb(dataTranspose, VbStructure, searchControl, trainControl, nthread, VB, HMM, HMMVB, bprint));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HDclust_rcpp_findModes", (DL_FUNC) &_HDclust_rcpp_findModes, 3},
    {"_HDclust_rcpp_clust", (DL_FUNC) &_HDclust_rcpp_clust, 5},
    {"_HDclust_rcpp_trainHmmVb", (DL_FUNC) &_HDclust_rcpp_trainHmmVb, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_HDclust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
