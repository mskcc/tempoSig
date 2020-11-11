#' NMF Iteration
#' 
#' Non-negative matrix factorization with either fixed or variable signature matrix.
#' 
#' If \code{denovo = TRUE}, the full NMF is performed with both \code{W} and \code{H} matrix 
#' determined de novo. If \code{denovo = FALSE}, \code{W} is fixed as input and only \code{H} 
#' is determined.
#' 
#' @param object \code{tempoSig}.
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
nmf <- function(catalog, denovo = FALSE, signat = NULL, K = NULL, nrun = 10, ard = FALSE,
                   verbose = FALSE, progress.bar = FALSE, Itmax = 100000, Tol = 1e-4,
                   a = 10, nprint = 100, penalizer='l2'){
  
  mat <- catalog
  if(!is.matrix(mat)) mat <- as.matrix(catalog)
  nullr <- sum(Matrix::rowSums(mat)==0)
  nullc <- sum(Matrix::colSums(mat)==0)
  if(nullr>0 & denovo) stop('Input matrix contains empty rows')
  if(nullc>0) stop('Input matrix contains empty columns')
  
  if(!penalizer %in% c('l1','l2')) stop('Unknown penalizer')
  
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  
  if(is.null(K)){
    if(!is.null(signat))
      K <- NCOL(signat)
    else
      K <- nrow
  }
  rank <- K
  
  if(ard){
    muv <- mean(mat)
    if(penalizer=='l1'){
      cc = nrow + ncol + a + 1  # ARD hyperparameters
      b = sqrt((a - 1)*(a - 2)*muv / K)
    }
    else{
      cc = (nrow + ncol)/2 + a + 1
      b = pi*(a - 1)*muv / (2*K)
    }
  }

  if(verbose) progress.bar <- FALSE  # can't show messages within progress.bar
  
  if(progress.bar) pb <- txtProgressBar(style = 3)
  rmax <- -Inf
  for(irun in seq(nrun)){
    if(verbose) cat('Run #',irun,':\n')
    wh <- init(nrow, ncol, mat, rank)
    if(!is.null(signat)){
      wh$ew <- signat
      rownames(wh$eh) <- colnames(wh$ew)
    }
    lkold <- -Inf
    tol <- Inf
    for(it in seq_len(Itmax)){
      if(ard){ 
        if(it > 1) lambda0 <- lambda
        if(penalizer=='l1')
          lambda <- (colSums(wh$ew) + rowSums(wh$eh) + b) / cc
        else
          lambda <- ((colSums(wh$ew^2) + rowSums(wh$eh)^2) / 2.0 + b) / cc
      } else lambda <- NULL
      wh <- nmf_updateR(mat, wh$ew, wh$eh, nrow, ncol, rank, denovo = denovo, ard = ard, 
                        lambda = lambda, penalizer = penalizer)
      lk0 <- likelihood(mat, wh$ew, wh$eh)
      if(ard){
        if(penalizer=='l1')
          lk0 <- lk0 - sum((colSums(wh$ew) + rowSums(wh$eh) + b) / lambda + cc*log(lambda))
        else
          lk0 <- lk0 - sum(((colSums(wh$ew^2) + rowSums(wh$eh^2)) / 2.0 + b) / lambda + cc*log(lambda))
      }
      if(it > 1){
        if(ard){ 
          tol <- max(abs((lambda - lambda0)/lambda0))
          if(tol < Tol) break
        } else if(abs(lkold - lk0) < Tol * abs(lkold)) break
      }
      if(verbose) if(it %% nprint == 0){
        if(ard){
          B <- b / cc
          Keff <- which((lambda - B)/B > tol)
        }
        cat(it,': likelihood = ',lk0, sep='')
        if(ard) cat(', K = ',length(Keff))
        cat('\n')
      }
      lkold <- lk0
    }
    if(ard){ 
      B <- b / cc
      Keff <- which((lambda - B)/B > tol)
    } else Keff <- K
    if(verbose){ 
      cat('Nsteps =', it, ', likelihood =', lk0, sep = '')
      if(ard) cat(', K = ', length(Keff), '\n', sep = '')
      else cat('\n',sep='')
    }
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
  
  if(ard){
    wmax <- wmax[, Kmax, drop = FALSE]
    hmax <- hmax[Kmax, , drop = FALSE]
  }
  cs <- colSums(wmax)
  wmax <- t(t(wmax) / cs)
  hmax <- hmax * cs
  signat <- wmax
  H <- t(t(hmax)/colSums(hmax))
  
  return(list(W = signat, H = H, logLik = rmax))
}

# single update step of NMF
nmf_updateR <- function(x, w, h, n, m, r, denovo = TRUE, ard = FALSE, lambda = NULL,
                        penalizer = 'l2'){
  
  x <- as.matrix(x)
  w <- as.matrix(w)
  h <- as.matrix(h)
  
  up <- h*(t(w) %*% (x/(w %*% h)))
  wi <- colSums(w)
  if(ard){ 
    if(penalizer=='l1')
      wi <- wi + 1.0 / lambda
    else
      wi <- wi + colSums(w) / lambda
  }
  down <- matrix(rep(wi, m), nrow = r, ncol = m, byrow = FALSE)
  
  h <- up / down
  h[h < .Machine$double.eps] <- .Machine$double.eps
  
  if(denovo){
    up <- w*((x / (w %*% h)) %*% t(h))
    hj <- rowSums(h)
    if(ard){ 
      if(penalizer=='l1')
        hj <- hj + 1.0 / lambda
      else
        hj <- hj + rowSums(h) / lambda
    }
    down <- matrix(rep(hj, n), nrow = n, ncol = r, byrow = TRUE)
    w <- up / down
    w[w < .Machine$double.eps] <- .Machine$double.eps
  }
  
  list(ew=w, eh=h)
}

# initialize w and h
init <- function(nrow, ncol, mat, rank, max = 1.0){
  w <- matrix(stats::runif(n = nrow * rank), nrow = nrow, ncol = rank)
  rownames(w) <- rownames(mat)
  colnames(w) <- seq_len(rank)
  h <- matrix(stats::runif(n = rank * ncol), nrow = rank, ncol = ncol)
  rownames(h) <- seq_len(rank)
  colnames(h) <- colnames(mat)
  list(ew = w, eh = h)
}

likelihood <- function(mat, w, h, per.element = FALSE){
  
  wh <- as.vector(w %*% h)
  amat <- as.vector(mat)
  x <- sum(amat * log(wh) - wh)
  z <- amat[amat > 0]
  x <- x + sum(-z * log(z) + z)
  if(per.element) x <- x / nrow(mat) / ncol(mat)
  
  return(x)
}