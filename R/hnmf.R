#' NMF Iteration
#' 
#' Non-negative matrix factorization with either fixed or variable signature matrix.
#' 
#' If \code{denovo = TRUE}, the full NMF is performed with both \code{W} and \code{H} matrix 
#' determined de novo. If \code{denovo = FALSE}, \code{W} is fixed as input and only \code{H} 
#' is determined.
#' 
#' @param catalog Catalog matrix
#' @param signat Signature matrix
#' @param denovo Full iteration of both signature and exposure; if \code{FALSE} (\code{hnmf}),
#'        only the exposure is fit (requires \code{signat} input).
#' @param K Number of signatures. 
#' @param nrun Number of independent runs to generate
#' @param verbose Verbosity level
#' @param progress.bar Display progress bar.
#' @param Itmax Maximum no. of iteration. 
#' @param Tol Tolerance for checking convergence.
#'        
#' @return Object of class \code{tempoSig}.
#' 
#' @export
nmf <- function(catalog, sig = NULL, denovo = TRUE, K = NULL, nrun = 10, verbose = TRUE, 
                progress.bar = FALSE, Itmax = 100000, Tol = 1e-4, a = 10, nprint = 100,
                useC = FALSE){
  
  mat <- catalog
  if(!is.matrix(mat)) mat <- as.matrix(catalog)
  nullr <- sum(Matrix::rowSums(mat)==0)
  nullc <- sum(Matrix::colSums(mat)==0)
# if(nullr>0 & denovo) stop('Input matrix contains empty rows')
  if(nullc>0) stop('Input matrix contains empty columns')
  
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  
  if(is.null(K)) if(!is.null(sig)) K <- NCOL(sig)
  rank <- K
  fudge <- .Machine$double.eps  # to avoid zero denominator
#  fudge <- 0
     
  if(verbose) progress.bar <- FALSE  # can't show messages within progress.bar
  
  if(progress.bar) pb <- txtProgressBar(style = 3)
  rmax <- -Inf
  for(irun in seq(nrun)){
    if(verbose) cat('Run #',irun,':\n')
#   wh <- init(nrow, ncol, mat, rank)
    hyper <- list(aw = 0.0, ah = 0.0, bw = 1.0, bh = 1.0) # ML prior
    wh <- vb_init(nrow = nrow, ncol = ncol, mat = mat, rank = rank, hyper = hyper,
                  initializer = 'random')
    if(!is.null(sig)){
      wh$ew <- sig
      rownames(wh$eh) <- colnames(wh$ew)
    }
    lkold <- -Inf
    tol <- Inf
    for(it in seq_len(Itmax)){
      if(useC){
        wh <- vbnmf_update(mat, wh, hyper, c(fudge))
        wh$lw <- wh$w
        wh$lh <- wh$h
        lk0 <- likelihood(mat, wh$ew, wh$eh)
      }
      else{
        wh <- nmf_updateR(mat, wh$ew, wh$eh, nrow, ncol, rank, denovo = denovo)
        lk0 <- likelihood(mat, wh$ew, wh$eh)
      }
      if(it > 1) if(abs(lkold - lk0) < Tol * abs(lkold)) break
      if(verbose) if(it %% nprint == 0)
        cat(it,': likelihood = ', lk0, '\n', sep='')
      lkold <- lk0
    }
    Keff <- K
    if(verbose)
      cat('Nsteps =', it, ', likelihood =', lk0, '\n', sep = '')
    if((irun == 1 | lk0 > rmax) & !is.na(lk0)){
      rmax <- lk0
      wmax <- wh$ew
      hmax <- wh$eh
      Kmax <- Keff
    }
    if(verbose) cat('Max(likelihood) =',rmax,'\n')
    if(progress.bar) setTxtProgressBar(pb, irun/nrun)
  }
  if(progress.bar) close(pb)
  
  cs <- colSums(wmax)
  wmax <- t(t(wmax) / cs)
  hmax <- hmax * cs
  sig <- wmax
  H <- t(t(hmax)/colSums(hmax))
  
  return(list(W = sig, H = H, logLik = rmax))
}

# single update step of NMF
nmf_updateR <- function(x, w, h, n, m, r, denovo = TRUE){
  
  x <- as.matrix(x)
  w <- as.matrix(w)
  h <- as.matrix(h)
  
  if(denovo){
    up <- w*((x / (w %*% h)) %*% t(h))
    hj <- rowSums(h)
    down <- matrix(rep(hj, n), nrow = n, ncol = r, byrow = TRUE)
    w2 <- up / down
    w2[w2 < .Machine$double.eps] <- .Machine$double.eps
  }
  
  up <- h*(t(w) %*% (x/(w %*% h)))
  if(denovo) w <- w2
  
  wi <- colSums(w)
  down <- matrix(rep(wi, m), nrow = r, ncol = m, byrow = FALSE)
  h <- up / down
  h[h < .Machine$double.eps] <- .Machine$double.eps

  return(list(ew=w, eh=h))
}

# initialize w and h
init <- function(nrow, ncol, mat, rank, max = 1.0){
  w <- matrix(stats::runif(n = nrow * rank), nrow = nrow, ncol = rank)
  h <- matrix(stats::runif(n = rank * ncol), nrow = rank, ncol = ncol)
  rownames(w) <- rownames(mat)
  if(is.null(colnames(w))) 
    colnames(w) <- rownames(h) <- paste0('S',seq_len(rank))
  colnames(h) <- colnames(mat)
  
  return(list(ew = w, eh = h))
}

likelihood <- function(mat, w, h, per.element = FALSE){
  
  wh <- as.vector(w %*% h)
  amat <- as.vector(mat)
#  x <- sum(amat * log(wh) - wh)
  wh2 <- wh[wh > 0]
  x <- sum(amat[wh > 0] * log(wh2) -wh2)
  z <- amat[amat > 0]
  x <- x + sum(-z * log(z) + z)
  if(per.element) x <- x / nrow(mat) / ncol(mat)
  
  return(x)
}