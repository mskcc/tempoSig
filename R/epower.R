#' Estimate Power for Signature Detection
#' 
#' @param proportion Vector of signature proportions
#' @param reference Reference signature version; currently only supports \code{c('SA')}.
#' @param cut Remove signatures with proportions below \code{cut}.
#' @param use.model Use parametrized model. If \code{FALSE}, multinomial model will be simulated.
#' @param nsamp No. of samples to use for simulation
#' @param Msamp Vector mutation loads to sample
#' @param verbose Print messages during sampling
#' @param alpha False positive rate
#' @param Pmin Lower bound for power 
#' @export
#' 
powerEstimate <- function(proportion , reference = 'v2', use.model = TRUE, nsamp = 100,
                          Mcoeff = NULL, afit = NULL, bfit = NULL,
                          Msamp = NULL, verbose = TRUE, alpha = 0.05, Pmin = 0.8, ...){
  
  proportion <- proportion[proportion > 0]
  if(sum(proportion)!=1.0) proportion <- proportion / sum(proportion) # normalize if not
  
# Read reference signatures
  if(!reference %in% c('v2','SA')) stop('Reference other than v2 or SA not implemented')
  else if(reference=='SA')
    ref <- read.table(system.file('extdata', 'cosmic_SigAnalyzer_SBS_signatures.txt', 
                                 package = 'tempoSig'), header = TRUE, sep = '\t')
  else
    ref <- read.table(system.file('extdata', 'cosmic_snv_signatures_v2.txt', package ='tempoSig'),
                      header = TRUE, sep = '\t')
  
  sbs <- colnames(ref)
  if(reference=='v2') sbs <- colnames(ref) <- gsub('nature','',sbs)
  nsbs <- length(sbs)
  if(sum(is.na(names(proportion))) > 0){
    if(length(proportion) != nsbs) stop('Input proportion inconsistent with reference')
    names(proportion) <- sbs
  }
  else if(sum(!(names(proportion) %in% sbs)) > 0)
    stop('Input proportion inconsistent with reference')

  if(is.null(Mcoeff)){
    if(reference=='v2')
      mf <- read.table(system.file('extdata', 'Mcoeff_v2.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
    else
      mf <- read.table(system.file('extdata', 'Mcoeff.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
  } else mf <- Mcoeff
  
  etiol <- mf[,3]
  names(etiol) <- mf[,1]
    
  if(use.model){
  # Read parameter sets
    if(reference=='v2'){
      if(is.null(afit)) afit <- read.table(system.file('extdata', 'Acoeff_v2.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
      if(is.null(bfit)) bfit <- read.table(system.file('extdata', 'Bcoeff_v2.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
    } else{
      if(is.null(afit)) afit <- read.table(system.file('extdata', 'Acoeff.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
      if(is.null(bfit)) bfit <- read.table(system.file('extdata', 'Bcoeff.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
#    fmin <- read.table(system.file('extdata', 'Fmin.txt', package = 'tempoSig'), 
#                       header = TRUE, sep = '\t')
    }
    nsig <- length(sbs)
    if(nrow(afit)!=nsig | ncol(afit) != nsig | nrow(bfit)!=nsig | ncol(bfit)!=nsig)
      stop('afit/bfit dimension does not match signatures list')
    
    mf1 <- mf[,2]
    names(mf1) <- mf[,1]
  
    pr <- proportion   # predicted Mp
    for(i in seq_along(proportion)){
      sbs <- names(pr)[i]
      fb <- 1 - proportion[i]
      pb <- proportion[-i]/fb  # background proportion
      if(length(pb)==0)
        mpr <- mf1[sbs]
      else{
        a <- sum(afit[sbs, names(pb)]*pb)
        b <- sum(bfit[sbs, names(pb)]*pb)
#       fm <- max(fmin[sbs, names(pb)])
#       fm <- exp(sum(pb*log(fmin[sbs, names(pb)])))
        frq <- proportion[i]
#       if(frq < fm) mpr <- Inf
#       else{
        lp <- log10(frq)
        logm <- log10(mf1[sbs]) - a*lp + b*lp^2
        mpr <- 10^logm
      }
#     }
      pr[i] <- max(1,round(mpr))
    }
    dat <- data.frame(Proportion = proportion, Mp = pr, Etiology = etiol[names(pr)])
  } else{   # simulation
    if(is.null(Msamp)) Msamp <- c(outer(c(2, 5, 10), 10^seq(0,6), '*'))
    nt <- rownames(ref)
    pr <- rep(0, length(nt))
    names(pr) <- nt
    for(i in seq_along(proportion)){
      sbi <- names(proportion)[i]
      pr <- pr + ref[, sbi]*proportion[i] # mean catalog vector of signatures
    }
    
    x0 <- rep(0, nsbs)
    names(x0) <- sbs
    x0[names(proportion)[proportion > 0]] <- 1    # indicator for true positive
    
    for(iM in seq_along(Msamp)){
      M <- Msamp[iM]
      if(verbose) cat('M = ',M, '\n', sep='')
      xcat <- rmultinom(n = nsamp, size = M, prob = pr)
      rownames(xcat) <- nt
      colnames(xcat) <- seq(nsamp)
      ts <- tempoSig(data = xcat, signat=ref)
      ts <- extractSig(ts, compute.pval = TRUE, progress.bar = verbose, ...)
      ex <- expos(ts)
      pv <- pvalue(ts)
      tp <- rep(0, length(proportion))
      names(tp) <- names(proportion)
      for(i in seq(nsamp)){
        ss <- senspec(xhat = ex[i,], x = x0, pvalue=pv[i, ], alpha=alpha)
        tp[ss$tp] <- tp[ss$tp] + 1
      }
      tp <- as.data.frame(t(tp) / nsamp)
      pw <- cbind(data.frame(M=M), tp)
      if(iM==1) dpw <- pw
      else dpw <- rbind(dpw, pw)
    }
    pr <- proportion   # predicted Mp
    for(i in seq_along(dpw[, -1])){
      pw <- dpw[, i+1]
      if(max(pw) < Pmin) mpr <- Inf
      else{
        df <- min(10, length(pw))
        sp <- smooth.spline(x=log10(dpw[, 1]), y=pw, df=df)
        xpr <- predict(sp, x=seq(log10(dpw[1, 1]), log10(dpw[NROW(dpw), 1]), 0.01))
        idx <- which(min(abs(xpr$y - Pmin))==abs(xpr$y - Pmin))
        mpr <- 10^xpr$x[idx]
      }
      pr[i] <- round(mpr)
    }
    dat <- list(Mpr = data.frame(Proportion = proportion, Mp = pr, 
                                 Etiology = etiol[names(pr)]), pcurve = dpw)
  }
  return(dat)
}
