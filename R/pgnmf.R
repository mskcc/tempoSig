#' piggyback NMF
#' 
#' Enables de novo inference of single sample using pre-inferred exposure 
#' and signature sets for a large data set (stock)
#' @param x Catalog matrix to be analyzed
#' @param stock \code{tempoSig} object containing the stock solution
#' @export
pgnmf <- function(x, stock, progress.bar = TRUE){
  
  stock.x <- catalog(stock)
  stock.h <- expos(stock)
  stock.w <- signat(stock)
  sig <- rownames(stock.h)
  K <- nrow(stock.h)
  N <- ncol(stock.x)
  
  m <- colSums(x)
  n <- ncol(x)
  
#  dat <- do.call(rbind, lapply(seq(n), function(k){
  dat <- NULL
  if(progress.bar) pb <- txtProgressBar(style = 3)
  for(k in seq(n)){
    xcat <- cbind(stock.x, data.frame(xk = x[,k]))
    b <- tempoSig(data = xcat, signat = stock.w)
    expos(b) <- cbind(stock.h, data.frame(xk = rep(1/K, K)))
    tmb(b) <- c(tmb(b), 'xk' = m[k])
    b <- extractSig(b, method = 'bnmf', nrun = 1, Kmin = K, Kmax = K, initializer = 'restart',
                    verbose = 0)
    dx <- t(expos(b)[,N+1])
    rownames(dx) <- colnames(x)[k]
    dat <- rbind(dat, dx)
    if(progress.bar) setTxtProgressBar(pb, value = k/n)
#  }))
  }
  if(progress.bar) close(pb)
  return(dat)
}