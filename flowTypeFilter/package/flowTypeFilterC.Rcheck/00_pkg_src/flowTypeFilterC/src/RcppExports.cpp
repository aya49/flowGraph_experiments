// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// countCells
SEXP countCells(SEXP RFilter, SEXP RPartitionsPerMarker, SEXP RThresholds, SEXP RmaxPopSize, SEXP RthresChannels, SEXP RMfiData, SEXP RFrameExprData, SEXP RNumPops, SEXP Rverbose);
RcppExport SEXP _flowTypeFilterC_countCells(SEXP RFilterSEXP, SEXP RPartitionsPerMarkerSEXP, SEXP RThresholdsSEXP, SEXP RmaxPopSizeSEXP, SEXP RthresChannelsSEXP, SEXP RMfiDataSEXP, SEXP RFrameExprDataSEXP, SEXP RNumPopsSEXP, SEXP RverboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type RFilter(RFilterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RPartitionsPerMarker(RPartitionsPerMarkerSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RThresholds(RThresholdsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RmaxPopSize(RmaxPopSizeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RthresChannels(RthresChannelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RMfiData(RMfiDataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFrameExprData(RFrameExprDataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumPops(RNumPopsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rverbose(RverboseSEXP);
    rcpp_result_gen = Rcpp::wrap(countCells(RFilter, RPartitionsPerMarker, RThresholds, RmaxPopSize, RthresChannels, RMfiData, RFrameExprData, RNumPops, Rverbose));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP countCells(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_flowTypeFilterC_countCells", (DL_FUNC) &_flowTypeFilterC_countCells, 9},
    {"countCells",                  (DL_FUNC) &countCells,                  9},
    {NULL, NULL, 0}
};

RcppExport void R_init_flowTypeFilterC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
