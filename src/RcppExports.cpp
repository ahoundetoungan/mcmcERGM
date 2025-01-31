// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// updateOmega
arma::mat updateOmega(const double& dfa, const arma::mat& Sca, const int& M, const arma::mat& hete);
RcppExport SEXP _mcmcERGM_updateOmega(SEXP dfaSEXP, SEXP ScaSEXP, SEXP MSEXP, SEXP heteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type dfa(dfaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sca(ScaSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type hete(heteSEXP);
    rcpp_result_gen = Rcpp::wrap(updateOmega(dfa, Sca, M, hete));
    return rcpp_result_gen;
END_RCPP
}
// propdnorm
double propdnorm(const Eigen::VectorXd& x, const Eigen::VectorXd& mu, const Eigen::MatrixXd& invV);
RcppExport SEXP _mcmcERGM_propdnorm(SEXP xSEXP, SEXP muSEXP, SEXP invVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type invV(invVSEXP);
    rcpp_result_gen = Rcpp::wrap(propdnorm(x, mu, invV));
    return rcpp_result_gen;
END_RCPP
}
// propdnorm_eachm
Eigen::RowVectorXd propdnorm_eachm(const Eigen::MatrixXd& x, const Eigen::MatrixXd& V);
RcppExport SEXP _mcmcERGM_propdnorm_eachm(SEXP xSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(propdnorm_eachm(x, V));
    return rcpp_result_gen;
END_RCPP
}
// propdproposal
Eigen::RowVectorXd propdproposal(const Eigen::ArrayXXd& x, const Eigen::MatrixXd& V, const Eigen::RowVectorXd& js);
RcppExport SEXP _mcmcERGM_propdproposal(SEXP xSEXP, SEXP VSEXP, SEXP jsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::RowVectorXd& >::type js(jsSEXP);
    rcpp_result_gen = Rcpp::wrap(propdproposal(x, V, js));
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
// fsimhete
arma::mat fsimhete(const arma::mat& mu, const arma::mat& Sigma, const arma::rowvec& js, const int& khete, const int& M);
RcppExport SEXP _mcmcERGM_fsimhete(SEXP muSEXP, SEXP SigmaSEXP, SEXP jsSEXP, SEXP kheteSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type js(jsSEXP);
    Rcpp::traits::input_parameter< const int& >::type khete(kheteSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(fsimhete(mu, Sigma, js, khete, M));
    return rcpp_result_gen;
END_RCPP
}
// ffindcom
Eigen::MatrixXd ffindcom(const int& n);
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
arma::mat futility(const arma::cube& X, const arma::vec& theta, const double& hetval, const int& npar, const int& n, const bool& intercept);
RcppExport SEXP _mcmcERGM_futility(SEXP XSEXP, SEXP thetaSEXP, SEXP hetvalSEXP, SEXP nparSEXP, SEXP nSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type hetval(hetvalSEXP);
    Rcpp::traits::input_parameter< const int& >::type npar(nparSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(futility(X, theta, hetval, npar, n, intercept));
    return rcpp_result_gen;
END_RCPP
}
// fQrsym
double fQrsym(const Eigen::ArrayXXd& ar, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& wr, const int& npu, const int& npw, const int& nr);
RcppExport SEXP _mcmcERGM_fQrsym(SEXP arSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(fQrsym(ar, ur, wr, npu, npw, nr));
    return rcpp_result_gen;
END_RCPP
}
// fQrdir
double fQrdir(const Eigen::ArrayXXd& ar, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& vr, const Eigen::ArrayXXd& wr, const int& npu, const int& npv, const int& npw, const int& nr);
RcppExport SEXP _mcmcERGM_fQrdir(SEXP arSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npv(npvSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(fQrdir(ar, ur, vr, wr, npu, npv, npw, nr));
    return rcpp_result_gen;
END_RCPP
}
// fGibbsym
Eigen::ArrayXXd fGibbsym(const Eigen::ArrayXXd& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& wr, const int& npu, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbsym(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbsym(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fGibbdir
Eigen::ArrayXXd fGibbdir(const Eigen::ArrayXXd& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& vr, const Eigen::ArrayXXd& wr, const int& npu, const int& npv, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbdir(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
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
List fGibbsym2(const Eigen::ArrayXXd& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& wr, const int& npu, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbsym2(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
    Rcpp::traits::input_parameter< const int& >::type npu(npuSEXP);
    Rcpp::traits::input_parameter< const int& >::type npw(npwSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fGibbsym2(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R));
    return rcpp_result_gen;
END_RCPP
}
// fGibbdir2
List fGibbdir2(const Eigen::ArrayXXd& ar, const int& nblock, const int& ncombr, const arma::mat& combr, const arma::uvec& idrows, const arma::uvec& idcols, const arma::uvec& ident, const arma::uvec& ztncombr, const Eigen::ArrayXXd& ur, const Eigen::ArrayXXd& vr, const Eigen::ArrayXXd& wr, const int& npu, const int& npv, const int& npw, const int& nr, const int& R);
RcppExport SEXP _mcmcERGM_fGibbdir2(SEXP arSEXP, SEXP nblockSEXP, SEXP ncombrSEXP, SEXP combrSEXP, SEXP idrowsSEXP, SEXP idcolsSEXP, SEXP identSEXP, SEXP ztncombrSEXP, SEXP urSEXP, SEXP vrSEXP, SEXP wrSEXP, SEXP npuSEXP, SEXP npvSEXP, SEXP npwSEXP, SEXP nrSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ar(arSEXP);
    Rcpp::traits::input_parameter< const int& >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncombr(ncombrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combr(combrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idrows(idrowsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idcols(idcolsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ident(identSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ztncombr(ztncombrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type ur(urSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type vr(vrSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type wr(wrSEXP);
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
// fupdate_jstheta
double fupdate_jstheta(const double& jscal, const double& accept, const int& iteration, const double& target, const double& kappa, const double& jmin, const double& jmax);
RcppExport SEXP _mcmcERGM_fupdate_jstheta(SEXP jscalSEXP, SEXP acceptSEXP, SEXP iterationSEXP, SEXP targetSEXP, SEXP kappaSEXP, SEXP jminSEXP, SEXP jmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type jscal(jscalSEXP);
    Rcpp::traits::input_parameter< const double& >::type accept(acceptSEXP);
    Rcpp::traits::input_parameter< const int& >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< const double& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double& >::type jmin(jminSEXP);
    Rcpp::traits::input_parameter< const double& >::type jmax(jmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(fupdate_jstheta(jscal, accept, iteration, target, kappa, jmin, jmax));
    return rcpp_result_gen;
END_RCPP
}
// fupdate_jshete
arma::rowvec fupdate_jshete(const arma::rowvec& jscal, const arma::rowvec& accept, const int& iteration, const double& target, const double& kappa, const double& jmin, const double& jmax);
RcppExport SEXP _mcmcERGM_fupdate_jshete(SEXP jscalSEXP, SEXP acceptSEXP, SEXP iterationSEXP, SEXP targetSEXP, SEXP kappaSEXP, SEXP jminSEXP, SEXP jmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type jscal(jscalSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type accept(acceptSEXP);
    Rcpp::traits::input_parameter< const int& >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< const double& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double& >::type jmin(jminSEXP);
    Rcpp::traits::input_parameter< const double& >::type jmax(jmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(fupdate_jshete(jscal, accept, iteration, target, kappa, jmin, jmax));
    return rcpp_result_gen;
END_RCPP
}
// frecentering
arma::mat frecentering(arma::vec& theta, const arma::mat& hetval, const arma::uvec intindex);
RcppExport SEXP _mcmcERGM_frecentering(SEXP thetaSEXP, SEXP hetvalSEXP, SEXP intindexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type hetval(hetvalSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type intindex(intindexSEXP);
    rcpp_result_gen = Rcpp::wrap(frecentering(theta, hetval, intindex));
    return rcpp_result_gen;
END_RCPP
}
// frVtoM
List frVtoM(const Eigen::VectorXd& u, const Rcpp::IntegerVector& N, const double& M);
RcppExport SEXP _mcmcERGM_frVtoM(SEXP uSEXP, SEXP NSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(frVtoM(u, N, M));
    return rcpp_result_gen;
END_RCPP
}
// frVtoMnorm
List frVtoMnorm(const arma::vec& u, const IntegerVector& N, const double& M);
RcppExport SEXP _mcmcERGM_frVtoMnorm(SEXP uSEXP, SEXP NSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(frVtoMnorm(u, N, M));
    return rcpp_result_gen;
END_RCPP
}
// frMtoV
Eigen::VectorXd frMtoV(List& u, const Rcpp::IntegerVector& N, const double& M);
RcppExport SEXP _mcmcERGM_frMtoV(SEXP uSEXP, SEXP NSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(frMtoV(u, N, M));
    return rcpp_result_gen;
END_RCPP
}
// frMceiltoV
Eigen::VectorXd frMceiltoV(List& u, const Rcpp::IntegerVector& N, const double& M);
RcppExport SEXP _mcmcERGM_frMceiltoV(SEXP uSEXP, SEXP NSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(frMceiltoV(u, N, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mcmcERGM_updateOmega", (DL_FUNC) &_mcmcERGM_updateOmega, 4},
    {"_mcmcERGM_propdnorm", (DL_FUNC) &_mcmcERGM_propdnorm, 3},
    {"_mcmcERGM_propdnorm_eachm", (DL_FUNC) &_mcmcERGM_propdnorm_eachm, 2},
    {"_mcmcERGM_propdproposal", (DL_FUNC) &_mcmcERGM_propdproposal, 3},
    {"_mcmcERGM_fsimtheta", (DL_FUNC) &_mcmcERGM_fsimtheta, 4},
    {"_mcmcERGM_fsimhete", (DL_FUNC) &_mcmcERGM_fsimhete, 5},
    {"_mcmcERGM_ffindcom", (DL_FUNC) &_mcmcERGM_ffindcom, 1},
    {"_mcmcERGM_fdatar", (DL_FUNC) &_mcmcERGM_fdatar, 4},
    {"_mcmcERGM_futility", (DL_FUNC) &_mcmcERGM_futility, 6},
    {"_mcmcERGM_fQrsym", (DL_FUNC) &_mcmcERGM_fQrsym, 6},
    {"_mcmcERGM_fQrdir", (DL_FUNC) &_mcmcERGM_fQrdir, 8},
    {"_mcmcERGM_fGibbsym", (DL_FUNC) &_mcmcERGM_fGibbsym, 14},
    {"_mcmcERGM_fGibbdir", (DL_FUNC) &_mcmcERGM_fGibbdir, 16},
    {"_mcmcERGM_fGibbsym2", (DL_FUNC) &_mcmcERGM_fGibbsym2, 14},
    {"_mcmcERGM_fGibbdir2", (DL_FUNC) &_mcmcERGM_fGibbdir2, 16},
    {"_mcmcERGM_fIDsym", (DL_FUNC) &_mcmcERGM_fIDsym, 2},
    {"_mcmcERGM_fIDdir", (DL_FUNC) &_mcmcERGM_fIDdir, 2},
    {"_mcmcERGM_fupdate_jstheta", (DL_FUNC) &_mcmcERGM_fupdate_jstheta, 7},
    {"_mcmcERGM_fupdate_jshete", (DL_FUNC) &_mcmcERGM_fupdate_jshete, 7},
    {"_mcmcERGM_frecentering", (DL_FUNC) &_mcmcERGM_frecentering, 3},
    {"_mcmcERGM_frVtoM", (DL_FUNC) &_mcmcERGM_frVtoM, 3},
    {"_mcmcERGM_frVtoMnorm", (DL_FUNC) &_mcmcERGM_frVtoMnorm, 3},
    {"_mcmcERGM_frMtoV", (DL_FUNC) &_mcmcERGM_frMtoV, 3},
    {"_mcmcERGM_frMceiltoV", (DL_FUNC) &_mcmcERGM_frMceiltoV, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mcmcERGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
