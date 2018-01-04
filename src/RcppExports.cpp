// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// bess_lm
List bess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, int T, int max_steps, Eigen::VectorXd beta0);
RcppExport SEXP _BeSS_bess_lm(SEXP XSEXP, SEXP ySEXP, SEXP TSEXP, SEXP max_stepsSEXP, SEXP beta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta0(beta0SEXP);
    rcpp_result_gen = Rcpp::wrap(bess_lm(X, y, T, max_steps, beta0));
    return rcpp_result_gen;
END_RCPP
}
// get_A
List get_A(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, double coef0, int T, Eigen::VectorXd B, Eigen::VectorXd weights);
RcppExport SEXP _BeSS_get_A(SEXP XSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP coef0SEXP, SEXP TSEXP, SEXP BSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_A(X, y, beta, coef0, T, B, weights));
    return rcpp_result_gen;
END_RCPP
}
// getcox_A
List getcox_A(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::VectorXd beta, int T, Eigen::VectorXd B, Eigen::VectorXd status, Eigen::VectorXd weights);
RcppExport SEXP _BeSS_getcox_A(SEXP XSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP TSEXP, SEXP BSEXP, SEXP statusSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type status(statusSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(getcox_A(X, y, beta, T, B, status, weights));
    return rcpp_result_gen;
END_RCPP
}
// EigenR
Eigen::MatrixXd EigenR(Eigen::MatrixXd X);
RcppExport SEXP _BeSS_EigenR(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR(X));
    return rcpp_result_gen;
END_RCPP
}
// gbess_lm
List gbess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd G, Eigen::VectorXd index, List PhiG, List invPhiG, int T0, int max_steps, Eigen::VectorXd beta0, int n, int p, int N);
RcppExport SEXP _BeSS_gbess_lm(SEXP XSEXP, SEXP ySEXP, SEXP GSEXP, SEXP indexSEXP, SEXP PhiGSEXP, SEXP invPhiGSEXP, SEXP T0SEXP, SEXP max_stepsSEXP, SEXP beta0SEXP, SEXP nSEXP, SEXP pSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type index(indexSEXP);
    Rcpp::traits::input_parameter< List >::type PhiG(PhiGSEXP);
    Rcpp::traits::input_parameter< List >::type invPhiG(invPhiGSEXP);
    Rcpp::traits::input_parameter< int >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(gbess_lm(X, y, G, index, PhiG, invPhiG, T0, max_steps, beta0, n, p, N));
    return rcpp_result_gen;
END_RCPP
}
// gget_A
List gget_A(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd G, Eigen::VectorXd index, int T0, Eigen::VectorXd beta0, double coef0, int n, int p, int N, Eigen::VectorXd weights, Eigen::VectorXd B00);
RcppExport SEXP _BeSS_gget_A(SEXP XSEXP, SEXP ySEXP, SEXP GSEXP, SEXP indexSEXP, SEXP T0SEXP, SEXP beta0SEXP, SEXP coef0SEXP, SEXP nSEXP, SEXP pSEXP, SEXP NSEXP, SEXP weightsSEXP, SEXP B00SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type index(indexSEXP);
    Rcpp::traits::input_parameter< int >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< double >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type B00(B00SEXP);
    rcpp_result_gen = Rcpp::wrap(gget_A(X, y, G, index, T0, beta0, coef0, n, p, N, weights, B00));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BeSS_bess_lm", (DL_FUNC) &_BeSS_bess_lm, 5},
    {"_BeSS_get_A", (DL_FUNC) &_BeSS_get_A, 7},
    {"_BeSS_getcox_A", (DL_FUNC) &_BeSS_getcox_A, 7},
    {"_BeSS_EigenR", (DL_FUNC) &_BeSS_EigenR, 1},
    {"_BeSS_gbess_lm", (DL_FUNC) &_BeSS_gbess_lm, 12},
    {"_BeSS_gget_A", (DL_FUNC) &_BeSS_gget_A, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_BeSS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
