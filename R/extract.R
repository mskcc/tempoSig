#' Infer Signature Proportions
#'
#' Use known signature list to find the most likely exposures in samples
#'
#' @param object Object of class \code{tempoSig}
#' @param method Refitting method; \code{mle} for maximum likelihood (default) or
#'               \code{mutCone} for mutationalCone.
#' @param itmax Maximum number of iterations for maximum likelihood estimate
#' @param tol Tolerance for convergence
#' @param min.tmb Minimum number of mutations in each sample. If \code{tmb} is less,
#'        \code{NA} will be returned for the sample.
#' @param compute.pval Estimate p-values
#' @param nperm Number of permutations
#' @param progress.bar Display progress bar
#' @import Rcpp
#' @useDynLib tempoSig
#' @examples
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt', package='tempoSig'))
#' b <- tempoSig(data)
#' b <- extractSig(b, progress.bar = TRUE)
#' b_pv <- extractSig(b, compute.pval = TRUE, progress.bar = TRUE)
#' @export
extractSig <- function(object, method = 'mle', itmax = 1000, tol = 1e-4, min.tmb = 5,
                       compute.pval = FALSE, nperm = 1000, progress.bar = FALSE){

  if(!is(object, 'tempoSig')) stop('object is not of class tempoSig')
  if(!method %in% c('mle','mutcone')) stop('Method is either mle or mutcone')

  spectrum <- catalog(object)
  ref <- signat(object)
  mut.load <- tmb(object)

  nt <- rownames(ref)
  nnt <- length(nt)
  idx <- match(nt, rownames(spectrum))
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Mutation types in spectrum do not match reference.')
  spectrum <- spectrum[idx, ]  # rearrange rows to match reference

  nsample <- length(tmb(object))
  nref <- NCOL(ref)
  h <- matrix(0, nrow = nsample, ncol = nref)
  colnames(h) <- colnames(ref)
  rownames(h) <- colnames(spectrum)
  if(compute.pval) pv <- h

  if(progress.bar) pb <- txtProgressBar(style = 3)
  for(i in seq(nsample)){
    spec <- spectrum[, i]
    if(method=='mle')
      h[i, ] <- hi <- fitMLE(x = spec, ref = ref, itmax = itmax, tol = tol, 
                             min.tmb = min.tmb)
    else
      h[i, ] <- hi <- mutationalCone(catalog = spectrum[, i, drop=F], 
                                     signature = ref, normalize = TRUE)
    if(compute.pval){
      perm <- matrix(1, nrow=nref, ncol=nperm)
      for(k in seq(nperm)){   # samples under null hypothesis
        rsp <- spec[sample(nnt)]
        names(rsp) <- nt
        if(method=='mle')
          perm[, k] <- fitMLE(x = rsp, ref = ref, itmax = itmax, tol = tol, min.tmb = min.tmb)
        else{
          rsp <- matrix(rsp, ncol=1)
          rownames(rsp) <- nt
          perm[, k] <- mutationalCone(catalog = rsp, signature = ref, normalize = TRUE)
        }
      }
      pv[i, ] <- rowSums(perm >= hi) / nperm
    }
    if(progress.bar) setTxtProgressBar(pb, i/nsample)
  }
  if(progress.bar) close(pb)

  expos(object) <- h
  if(compute.pval) pvalue(object) <- pv

  return(object)
}

#' Maximum likelihood inference of signature proportions
#' 
#' @param x Vector of observed proportions of mutation contexts
#' @param ref Matrix of reference signature proportions with mutation types in rows
#'            and signatures in columns
#' @param itmax Maximum number of iterations
#' @param tol Tolerance of convergence
#' @param min.tmb Minimum mutation counts; if lower than this, \code{NA} is returned
#' @return Vector of estimated signature exposure proportions
#' @export
fitMLE <- function(x, ref, itmax = 1000, tol = 1e-4, min.tmb = 5){

  num_sigs <- NCOL(ref)
  num_muts <- sum(x)
  if(num_muts < min.tmb){
    warning(paste0('Sample has less than ', min.tmb, ' mutations.'))
    return(rep(NA, num_muts))
  }
  x <- x / num_muts
  x0 <- gtools::rdirichlet(n = 1, alpha = rep(10, num_sigs)) #initial guess
  p <- mlestimate(x, x0, ref, Itmax=itmax, Tol=tol)
  h <- p$x^2/sum(p$x^2)
  names(h) <- colnames(ref)

  return(h)
}
