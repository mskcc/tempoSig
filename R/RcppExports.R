# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

mlestimate <- function(D, x0, ref, Itmax, Tol) {
    .Call('_tempoSig_mlestimate', PACKAGE = 'tempoSig', D, x0, ref, Itmax, Tol)
}

vbnmf_update <- function(X, wh, hyper, fudge) {
    .Call('_tempoSig_vbnmf_update', PACKAGE = 'tempoSig', X, wh, hyper, fudge)
}

