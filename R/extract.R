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
#' @param pvtest Algorithm for p-value computation; 
#'            \code{c('permutation','lrt','x.permutation')} for permutation resampling of 
#'            signatures, likelihood ratio test (asymptotic formula), or
#'            permutation of count data.
#' @param ... Other parameters for \code{denovo} with \code{method = 'hnmf'}
#' @import Rcpp
#' @useDynLib tempoSig
#' @examples
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt', package='tempoSig'))
#' b <- tempoSig(data)
#' b <- extractSig(b, progress.bar = TRUE)
#' b_pv <- extractSig(b, compute.pval = TRUE, progress.bar = TRUE)
#' @export
extractSig <- function(object, method = 'mle', itmax = 1000, tol = 1e-4, min.tmb = 2,
                       compute.pval = FALSE, nperm = 1000, progress.bar = FALSE,
                       pvtest = 'permutation', cosmic = FALSE, Kmin = 2, K = 2, 
                       Kmax = 30, nrun =10, useC = TRUE, initializer = 'random', ...){

  if(Kmax < 2 ) stop('Kmax must be at least 2')
  if(!is(object, 'tempoSig')) stop('object is not of class tempoSig')
  if(!method %in% c('nmf','hnmf','bnmf','mle','mutcone')) stop('Unknown method')
  if(method != 'mle' & pvtest == 'lrt') stop('Likelihood ratio test is possible only with MLE')
  if(method == 'mutcone' & compute.pval) stop('P-value in mutcone not implemented')

  spectrum <- catalog(object)
  ref <- signat(object)
  nref <- NCOL(ref)
  mut.load <- tmb(object)
  nsample <- length(tmb(object))
  
  nt <- rownames(ref)
  nnt <- length(nt)
  idx <- match(nt, rownames(spectrum))
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Mutation types in spectrum do not match reference.')
  spectrum <- spectrum[idx, , drop = FALSE]  # rearrange rows to match reference
  
  if(method %in% c('nmf','hnmf','bnmf')){
    if(compute.pval) cat('Computing exposures ...\n',sep='')
    if(method=='hnmf'){ 
      K <- NULL
      sigs <- ref
    } else if(method=='nmf') sigs <- NULL
    if(method %in% c('nmf','hnmf'))
      nmf <- nmf(catalog = spectrum, denovo = (method=='nmf'), sig = sigs, K = K,
               progress.bar = progress.bar, nrun = nrun, useC = useC, ...)
    else{ # bnmf
      sig <- colnames(signat(object))
      object <- bnmf(object, ranks=seq(Kmin,Kmax), nrun = nrun, useC = useC, 
                     initializer = initializer, ...)
      W <- signat(object)
      H <- expos(object)
      if(initializer=='restart')
        colnames(W) <- rownames(H) <- sig
      else
        colnames(W) <- rownames(H) <- paste0('S',seq(NCOL(W)))
      signat(object) <- W
      expos(object) <- H
      return(object)
    }
    
    H <- t(nmf$H)
    if(method != 'hnmf'){ 
      W <- nmf$W
      if(is.null(rownames(W))) rownames(W) <- nt
      if(is.null(colnames(W))) colnames(W) <- S <- paste0('S', seq(ncol(W)))
      if(is.null(colnames(H))) colnames(H) <- S
      if(is.null(rownames(H))) rownames(H) <- colnames(spectrum)
#     if(cosmic){
#       cosine <- cosineSimilarity(A = W, B = ref, diag = FALSE)
#       lsap <- clue::solve_LSAP(cosine, maximum = TRUE)
#       sbs <- colnames(cosine)[lsap]
#       colnames(W) <- colnames(H) <- sbs
#       misc(object) <- list(cosine = cosine)
#     }
      signat(object) <- W
    }
    expos(object) <- H
    logLik(object) <- nmf$logLik
    if(compute.pval & method == 'hnmf'){   # estimate p-values by permutation resampling
      cat('Estimating p-values ...\n',sep='')
      if(progress.bar) pb <- txtProgressBar(style = 3)
      pv <- matrix(0, nrow = nsample, ncol = nref)
      rownames(pv) <- colnames(spectrum)
      colnames(pv) <- colnames(ref)
      for(l in seq(nperm)){
        if(pvtest == 'x.permutation'){
          spec <- apply(spectrum, 2, function(x){sample(x, size = length(x), replace = FALSE)})
          rownames(spec) <- rownames(spectrum)
          nmf <- nmf(catalog = spec, denovo = FALSE, sig = ref, progress.bar = FALSE, 
                   nrun = 1, verbose = FALSE, ...)
          pv <- pv + (t(nmf$H) >= H)
        } else if(pvtest == 'permutation'){
          rref <- ref
          for(k in seq(nref)){
            rref[, k] <- rref[sample(nnt),k]
            names(rref[,k]) <- nt
            nmf <- nmf(catalog = spectrum, denovo = FALSE, sig = rref, progress.bar = FALSE, 
                       nrun = 1, verbose = FALSE, ...)
            pv[, k] <- pv[, k] + (t(nmf$H)[, k] >= H[, k])
          }
        }
        if(progress.bar) setTxtProgressBar(pb, l/nperm)
      }
      pvalue(object) <- pv/nperm
      if(progress.bar) close(pb)
    }
    return(object)
  }
  
  h <- matrix(0, nrow = nsample, ncol = nref)
  signatures <- colnames(ref)
  colnames(h) <- signatures
  rownames(h) <- colnames(spectrum)
  if(compute.pval) pv <- h
  
  if(progress.bar) pb <- txtProgressBar(style = 3)
  L <- rep(0, nsample)
  names(L) <- colnames(spectrum)
  for(i in seq(nsample)){
    if(mut.load[i] < min.tmb){ 
      h[i, ] <- rep(NA, nref)
      if(compute.pval) pv[i, ] <- rep(NA, nref)
      next()
    }
    spec <- spectrum[, i]
    if(method=='mle'){
      fi <- fitMLE(x = spec, ref = ref, itmax = itmax, tol = tol)
      h[i, ] <- hi <- fi$h
      L[i] <- fi$loglik
    }
    else
      h[i, ] <- hi <- mutationalCone(catalog = spectrum[, i, drop=F], 
                                     signature = ref, normalize = TRUE)
    if(compute.pval){
      if(pvtest == 'lrt'){
        for(k in seq(nref)){
          fm <- fitMLE(x = spec, ref = ref[,-k], itmax = itmax, tol = tol)
          q = 2*mut.load[i]*(L[i] - fm$loglik)
          pv[i, k] <- pchisq(q, df = 1, lower.tail = FALSE)
        }
      } else{
        perm <- matrix(0, nrow=nref, ncol=nperm)
        rownames(perm) <- signatures
        for(l in seq(nperm)){   # samples under null hypothesis
          rref <- ref
          if(pvtest == 'x.permutation'){
            rsp <- spec
            rsp <- rsp[sample(length(rsp))]
            names(rsp) <- nt
            fh <- fitMLE(x = rsp, ref = ref, itmax = itmax, tol = tol)
            perm[, l] <- fh$h
            if(progress.bar) if(l %% 100 == 0)
              setTxtProgressBar(pb, ((i - 1)*nperm + l)/nsample/nperm)
          } else if(pvtest == 'permutation'){
            for(k in seq(nref)){
              rref[, k] <- rref[sample(nnt),k]
              names(rref[,k]) <- nt
              fk <- fitMLE(x = spec, ref = rref, itmax = itmax, tol = tol)
              perm[k, l] <- fk$h[k]
              if(progress.bar) if(l %% 10 == 0)
                setTxtProgressBar(pb, ((i - 1)*nperm*nref + (l - 1)*nref + k)/nsample/nperm/nref)
            }
          } else stop('Unknown pvtest')
        }
        pv[i, ] <- rowSums(perm >= h[i, ]) / nperm
      }
    }
    if(progress.bar) setTxtProgressBar(pb, i/nsample)
  }
  if(progress.bar) close(pb)

  expos(object) <- h
  if(compute.pval) pvalue(object) <- pv
  if(method=='mle')
    logLik(object) <- L

  return(object)
}

#' Maximum likelihood inference of signature proportions
#' 
#' @param x Vector of observed proportions of mutation contexts
#' @param ref Matrix of reference signature proportions with mutation types in rows
#'            and signatures in columns
#' @param itmax Maximum number of iterations
#' @param tol Tolerance of convergence
#' @return List of exposure vector \code{p} and \code{loglik}
#' @export
fitMLE <- function(x, ref, itmax = 1000, tol = 1e-4){

  num_sigs <- NCOL(ref)
  num_muts <- sum(x)
  x <- x / num_muts
  x0 <- gtools::rdirichlet(n = 1, alpha = rep(10, num_sigs)) #initial guess
  p <- mlestimate(x, x0, ref, Itmax=itmax, Tol=tol)
  h <- p$x^2/sum(p$x^2)
  names(h) <- colnames(ref)
  
  return(list(h=h, loglik = p$lkh))
}
