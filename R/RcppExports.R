# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

propdnorm <- function(x, mu, invV) {
    .Call(`_mcmcERGM_propdnorm`, x, mu, invV)
}

fsimtheta <- function(mu, Sigma, js, npar) {
    .Call(`_mcmcERGM_fsimtheta`, mu, Sigma, js, npar)
}

ffindcom <- function(n) {
    .Call(`_mcmcERGM_ffindcom`, n)
}

fdatar <- function(X, ftovar, nvar, K) {
    .Call(`_mcmcERGM_fdatar`, X, ftovar, nvar, K)
}

futility <- function(X, theta, npar, n, intercept) {
    .Call(`_mcmcERGM_futility`, X, theta, npar, n, intercept)
}

fQrsym <- function(ar, ur, wr, npu, npw, nr) {
    .Call(`_mcmcERGM_fQrsym`, ar, ur, wr, npu, npw, nr)
}

fQrdir <- function(ar, ur, vr, wr, npu, npv, npw, nr) {
    .Call(`_mcmcERGM_fQrdir`, ar, ur, vr, wr, npu, npv, npw, nr)
}

fGibbsym <- function(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R) {
    .Call(`_mcmcERGM_fGibbsym`, ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R)
}

fGibbdir <- function(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R) {
    .Call(`_mcmcERGM_fGibbdir`, ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R)
}

fGibbsym2 <- function(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R) {
    .Call(`_mcmcERGM_fGibbsym2`, ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr, R)
}

fGibbdir2 <- function(ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R) {
    .Call(`_mcmcERGM_fGibbdir2`, ar, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr, R)
}

fIDsym <- function(M, nvec) {
    .Call(`_mcmcERGM_fIDsym`, M, nvec)
}

fIDdir <- function(M, nvec) {
    .Call(`_mcmcERGM_fIDdir`, M, nvec)
}

frMtoV <- function(u, N, M) {
    .Call(`_mcmcERGM_frMtoV`, u, N, M)
}

frMceiltoV <- function(u, N, M) {
    .Call(`_mcmcERGM_frMceiltoV`, u, N, M)
}

