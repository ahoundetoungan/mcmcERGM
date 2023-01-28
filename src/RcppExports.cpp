// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// propdnorm
double propdnorm(const arma::vec& x, const arma::vec& mu, const arma::mat& invV);
RcppExport SEXP _mcmcERGM_propdnorm(SEXP xSEXP, SEXP muSEXP, SEXP invVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invV(invVSEXP);
    rcpp_result_gen = Rcpp::wrap(propdnorm(x, mu, invV));
    return rcpp_result_gen;
END_RCPP
}
// fsimtheta
arma::vec fsimtheta(const arma::vec& mu, const arma::mat& Sigma, const double& js, const int& npar);
RcppExport SEXP _mcmcERGM_fsimtheta(SEXP muSEXP, SEXP SigmaSEXP, SEXP jsSEXP, SEXP nparSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type js(jsSEXP);
    Rcpp::traits::input_parameter< const int& >::type npar(nparSEXP);
    rcpp_result_gen = Rcpp::wrap(fsimtheta(mu, Sigma, js, npar));
    return rcpp_result_gen;
END_RCPP
}
// ffindcom
arma::mat ffindcom(const int& n);
RcppExport SEXP _mcmcERGM_ffindcom(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(ffindcom(n));
    return rcpp_result_gen;
END_RCPP
}
// fdatar
arma::cube fdatar(const arma::mat X, List ftovar, const int& nvar, const int& K);
RcppExport SEXP _mcmcERGM_fdatar(SEXP XSEXP, SEXP ftovarSEXP, SEXP nvarSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type ftovar(ftovarSEXP);
    Rcpp::traits::input_parameter< const int& >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(fdatar(X, ftovar, nvar, K));
    return rcpp_result_gen;
END_RCPP
}
// futility
arma::mat futility(const arma::cube& X, const arma::vec& theta, const int& npar, const int& n, const int& intercept);
RcppExport SEXP _mcmcERGM_futility(SEXP XSEXP, SEXP thetaSEXP, SEXP nparSEXP, SEXP nSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type npar(nparSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(futility(X, theta, npar, n, intercept));
    return rcpp_result_gen;
END_RCPP
}
// fQrsym
double fQrsym(const arma::mat& ar, const arma::mat& ur, const arma::mat& wr, const int& npu, const int& npw, const int& nr);
RcppExport SEXP _mcmcERGM_fQrsym(SEXP arSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(fQrsym(ar, ur, wr, npu, npw, nr));
    return rcpp_result_gen;
END_RCPP
}
// fQrdir
double fQrdir(const arma::mat& ar, const arma::mat& ur, const arma::mat& vr, const arma::mat& wr, const int& npu, const int& npv, const int& npw, const int& nr);
RcppExport SEXP _mcmcERGM_fQrdir(SEXP arSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npv(npvSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(fQrdir(ar, ur, vr, wr, npu, npv, npw, nr));
    return rcpp_result_gen;
END_RCPP
}
// fGibbsym
arma::mat fGibbsym(const arma::mat& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const arma::mat& ur, const arma::mat& wr, const int& npu, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbsym(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbsym(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fGibbdir
arma::mat fGibbdir(const arma::mat& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const arma::mat& ur, const arma::mat& vr, const arma::mat& wr, const int& npu, const int& npv, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbdir(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npv(npvSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbdir(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fGibbsym2
List fGibbsym2(const arma::mat& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const arma::mat& ur, const arma::mat& wr, const int& npu, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbsym2(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbsym2(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fGibbdir2
List fGibbdir2(const arma::mat& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const arma::mat& ur, const arma::mat& vr, const arma::mat& wr, const int& npu, const int& npv, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbdir2(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npv(npvSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbdir2(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fIDsym
List fIDsym(const int& M, const arma::vec& nvec);
RcppExport SEXP _mcmcERGM_fIDsym(SEXP MSEXP, SEXP nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nvec(nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(fIDsym(M, nvec));
    return rcpp_result_gen;
END_RCPP
}
// fIDdir
List fIDdir(const int& M, const arma::vec& nvec);
RcppExport SEXP _mcmcERGM_fIDdir(SEXP MSEXP, SEXP nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nvec(nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(fIDdir(M, nvec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mcmcERGM_propdnorm", (DL_FUNC) &_mcmcERGM_propdnorm, 3},
    {"_mcmcERGM_fsimtheta", (DL_FUNC) &_mcmcERGM_fsimtheta, 4},
    {"_mcmcERGM_ffindcom", (DL_FUNC) &_mcmcERGM_ffindcom, 1},
    {"_mcmcERGM_fdatar", (DL_FUNC) &_mcmcERGM_fdatar, 4},
    {"_mcmcERGM_futility", (DL_FUNC) &_mcmcERGM_futility, 5},
    {"_mcmcERGM_fQrsym", (DL_FUNC) &_mcmcERGM_fQrsym, 6},
    {"_mcmcERGM_fQrdir", (DL_FUNC) &_mcmcERGM_fQrdir, 8},
    {"_mcmcERGM_fGibbsym", (DL_FUNC) &_mcmcERGM_fGibbsym, 14},
    {"_mcmcERGM_fGibbdir", (DL_FUNC) &_mcmcERGM_fGibbdir, 16},
    {"_mcmcERGM_fGibbsym2", (DL_FUNC) &_mcmcERGM_fGibbsym2, 14},
    {"_mcmcERGM_fGibbdir2", (DL_FUNC) &_mcmcERGM_fGibbdir2, 16},
    {"_mcmcERGM_fIDsym", (DL_FUNC) &_mcmcERGM_fIDsym, 2},
    {"_mcmcERGM_fIDdir", (DL_FUNC) &_mcmcERGM_fIDdir, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mcmcERGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}