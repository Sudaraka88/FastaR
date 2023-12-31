// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extractAlnParam2
List extractAlnParam2(std::string file);
RcppExport SEXP _FastaR_extractAlnParam2(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(extractAlnParam2(file));
    return rcpp_result_gen;
END_RCPP
}
// extractAlnParam
List extractAlnParam(std::string file, int filter, double gap_thresh, double maf_thresh);
RcppExport SEXP _FastaR_extractAlnParam(SEXP fileSEXP, SEXP filterSEXP, SEXP gap_threshSEXP, SEXP maf_threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< double >::type gap_thresh(gap_threshSEXP);
    Rcpp::traits::input_parameter< double >::type maf_thresh(maf_threshSEXP);
    rcpp_result_gen = Rcpp::wrap(extractAlnParam(file, filter, gap_thresh, maf_thresh));
    return rcpp_result_gen;
END_RCPP
}
// getSNPs
std::vector<char> getSNPs(std::string file, int n_seq, int n_snp, std::vector<int> POS);
RcppExport SEXP _FastaR_getSNPs(SEXP fileSEXP, SEXP n_seqSEXP, SEXP n_snpSEXP, SEXP POSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n_seq(n_seqSEXP);
    Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type POS(POSSEXP);
    rcpp_result_gen = Rcpp::wrap(getSNPs(file, n_seq, n_snp, POS));
    return rcpp_result_gen;
END_RCPP
}
// extractSNPs
List extractSNPs(std::string file, int n_seq, int n_snp, std::vector<int> POS);
RcppExport SEXP _FastaR_extractSNPs(SEXP fileSEXP, SEXP n_seqSEXP, SEXP n_snpSEXP, SEXP POSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n_seq(n_seqSEXP);
    Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type POS(POSSEXP);
    rcpp_result_gen = Rcpp::wrap(extractSNPs(file, n_seq, n_snp, POS));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastaR_extractAlnParam2", (DL_FUNC) &_FastaR_extractAlnParam2, 1},
    {"_FastaR_extractAlnParam", (DL_FUNC) &_FastaR_extractAlnParam, 4},
    {"_FastaR_getSNPs", (DL_FUNC) &_FastaR_getSNPs, 4},
    {"_FastaR_extractSNPs", (DL_FUNC) &_FastaR_extractSNPs, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastaR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
