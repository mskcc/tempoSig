#' piggyback NMF
#' 
#' Enables de novo inference of single sample using pre-inferred exposure 
#' and signature sets for a large data set (stock)
#' @param x Catalog matrix to be analyzed
#' @param stock \code{tempoSig} object containing the stock solution
#' @param filter Filter exposure using CV cutoff \code{cutoff}
#' @param cutoff CV cutoff matrix for filtering
#' @return List of \code{mean} and coefficient of variation \code{cv} (sd/mean) of exposure matrices 
#' @export
pgnmf <- function(x, stock = NULL, filter = FALSE, cutoff = NULL, progress.bar = TRUE, verbose = 1, ...){
  
  if(is.null(stock)){
    fl <- system.file('extdata/pig_K11.rds', package = 'tempoSig')
    if(!file.exists(fl)) stop('Default stock object does not exist')
    stock <- readRDS(fl)
  }
  stock.x <- catalog(stock)
  stock.h <- expos(stock)
  stock.w <- signat(stock)
  sig <- rownames(stock.h)
  K <- nrow(stock.h)
  N <- ncol(stock.x)
  
  m <- colSums(x)
  n <- ncol(x)
  
#  dat <- do.call(rbind, lapply(seq(n), function(k){
  dat <- vdat <- NULL
  if(progress.bar) pb <- txtProgressBar(style = 3)
  for(k in seq(n)){
    xcat <- cbind(stock.x, data.frame(xk = x[,k]))
    b <- tempoSig(data = xcat, signat = stock.w)
    expos(b) <- cbind(stock.h, data.frame(xk = rep(1/K, K)))
    tmb(b) <- c(tmb(b), 'xk' = m[k])
    b <- extractSig(b, method = 'bnmf', nrun = 1, Kmin = K, Kmax = K, initializer = 'restart',
                    verbose = verbose -1 , ...)
    dx <- t(expos(b)[,N+1])
    rownames(dx) <- colnames(x)[k]
    dat <- rbind(dat, dx)
    dv <- t(misc(b)$dH[[1]][,N+1])
    rownames(dv) <- colnames(x)[k]
    vdat <- rbind(vdat, dv)
    if(progress.bar) setTxtProgressBar(pb, value = k/n)
#  }))
  }
  if(progress.bar) close(pb)
  
  if(filter){
#   cat('Filtering exposure...')
    if(is.null(cutoff)){
      fl <- system.file('extdata/pig_K11_cutoff.txt', package = 'tempoSig')
      if(!file.exists(fl)) stop('Default stock object does not exist')
      cutoff <- read.table(fl, header = TRUE, sep = '\t')
    }
    sig <- colnames(cutoff)[-1]
    if(!all(sig == colnames(dat)))
      stop('Signature names in cutoff do not agree with exposure matrix')
    sfit <- function(m, S){
      x <- log10(m)
      if(!S %in% sig) stop(paste0(S, ' not in signature set'))
      g <- smooth.spline(x = log10(cutoff$M), y = cutoff[, S], df = nrow(cutoff))
      gamma <- predict(g, x = x)$y
      return(gamma)
    }
    
    for(i in seq(n)) for(k in seq(K)){
      g <- sfit(m = m[i], S = sig[k])
      if(vdat[i, k] >= g) dat[i, k] <- 0
    }
  }
  
  return(list(mean=dat, cv=vdat))
}