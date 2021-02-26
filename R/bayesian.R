# Update hyperparameters
hyper_update <- function(hyper.update, wh, hyper, hyper.tying='ikj',
                         Niter=100, Tol=1e-4){
    
   if(sum(hyper.update)==0) return(hyper)
  
   aw0 <- hyper$aw
   ah0 <- hyper$ah
   if(hyper.tying=='ikj'){
     lwm <- mean(log(wh$lw))
     lhm <- mean(log(wh$lh))
     ewm <- mean(wh$ew)
     ehm <- mean(wh$eh)
   } else{
     lwm <- colMeans(log(wh$lw))
     lhm <- rowMeans(log(wh$lh))
     ewm <- colMeans(wh$ew)
     ehm <- rowMeans(wh$eh)
   }
   bw0 <- hyper$bw
   bh0 <- hyper$bh
  
   if(hyper.update[1]+hyper.update[3]>0){
    i <- 1
    while(i<Niter){
      if(hyper.update[1])
        dw <- (log(aw0)-digamma(aw0)-ewm/bw0+1+lwm-log(bw0))/
          (1/aw0-psigamma(aw0,1))
      else dw <- 0
      if(hyper.update[3])
        dh <- (log(ah0)-digamma(ah0)-ehm/bh0+1+lhm-log(bh0))/
          (1/ah0-psigamma(ah0,1))
      else dh <- 0
      aw1 <- aw0 - dw 
      ah1 <- ah0 - dh
      while(any(aw1<=0)){
        dw <- dw/2
        aw1 <- aw0 - dw 
      }
      while(any(ah1<=0)){
        dh <- dh/2
        ah1 <- ah0 - dh
      }
      
      df <- mean((1-aw1/aw0)^2)+mean((1-ah1/ah0)^2)
      if(df<Tol) break
      aw0 <- aw1
      ah0 <- ah1
      i <- i+1
    }
    if(i==Niter) stop('Hyper-parameter update failed to converge')
   } else{
     aw1 <- aw0
     ah1 <- ah0
   }
   if(hyper.update[2]) bw1 <- ewm
   else bw1 <- bw0
   if(hyper.update[4]) bh1 <- ehm
   else bh1 <- ehm
   list(aw=aw1, bw=bw1, ah=ah1, bh=bh1)
}

# Single update step in Bayesian NMF inference
vbnmf_updateR <- function(x, wh, r, hyper, fudge=NULL){

   x <- as.matrix(x)
   n <- dim(x)[1]
   m <- dim(x)[2]

   lw <- as.matrix(wh$lw)
   lh <- as.matrix(wh$lh)
   ew <- as.matrix(wh$ew)
   eh <- as.matrix(wh$eh)
   aw <- hyper$aw
   bw <- hyper$bw
   ah <- hyper$ah
   bh <- hyper$bh
  
   wth <- lw %*% lh
   sw <- lw*((x/wth) %*% t(lh))
   sh <- lh*(t(lw) %*% (x/wth))
  
#  alw <- aw + sw
#  bew <- 1/(aw/bw + t(replicate(n, rowSums(eh))))
   alw <- matrix(aw, nrow=n, ncol=r, byrow=TRUE) + sw
   bew <- 1/(matrix(aw/bw, nrow=n, ncol=r, byrow=TRUE) + 
                t(replicate(n, rowSums(eh))))
   ew <- alw*bew             # this update needs to precede lines below
  
#  alh <- ah + sh
#  beh <- 1/(ah/bh + replicate(m, colSums(ew)))
   alh <- matrix(ah, nrow=r, ncol=m, byrow=FALSE) + sh
   beh <- 1/(matrix(ah/bh, nrow=r, ncol=m, byrow=FALSE) + 
                replicate(m, colSums(ew)))
   eh <- alh*beh

   lw <- exp(digamma(alw))*bew
   lh <- exp(digamma(alh))*beh
   if(is.null(fudge)) fudge <- .Machine$double.eps
   lw[lw < fudge] <- fudge
   lh[lh < fudge] <- fudge
    
   wth <- lw %*% lh
   U1 <- -ew %*% eh - lgamma(x+1) - x*((((lw*log(lw))%*%lh) + 
          lw%*%(lh*log(lh)))/wth - log(wth))
   U2 <- -(aw/bw)*ew - lgamma(aw) + aw*log(aw/bw) + 
    alw*(1+log(bew))+lgamma(alw)
   U3 <- -(ah/bh)*eh - lgamma(ah) + ah*log(ah/bh) + 
    alh*(1+log(beh))+lgamma(alh)
   U <- sum(U1) + sum(U2) + sum(U3)
   U <- U/(n*m)  # log evidence per feature per cell
  
   w <- ew    
   h <- eh
   
   dw <- alw*bew^2
   dh <- alh*beh^2
    
   list(w=w, h=h, lw=lw, lh=lh, ew=ew, eh=eh, lkh=U, dw=dw, dh=dh)
}

# Initialize bNMF inference
vb_init <- function(nrow,ncol,mat,rank, max=1.0, hyper, initializer){
  
   if(initializer=='random'){
     w <- matrix(stats::rgamma(n=nrow*rank, shape=hyper$aw, 
           scale=hyper$bw/hyper$aw), nrow=nrow, ncol=rank, byrow=TRUE)
     h <- matrix(stats::rgamma(n=rank*ncol, shape=hyper$ah, 
           scale=hyper$bh/hyper$ah), nrow=rank, ncol=ncol, byrow=FALSE)
   }else stop('Unknown initializer')
  
   dw <- matrix(0, nrow=nrow, ncol=rank)
   dh <- matrix(0, nrow=rank, ncol=ncol)

   rownames(w) <- rownames(dw) <- rownames(mat)
   colnames(w) <- colnames(dw) <- seq_len(rank)
   rownames(h) <- rownames(dh) <- seq_len(rank)
   colnames(h) <- colnames(dh) <- colnames(mat)
  
   list(w=w, h=h, lw=w, lh=h, ew=w, eh=h, dw=dw, dh=dh)
}

#' Bayesian NMF inference of count matrix
#' 
#' Perform variational Bayes NMF and store factor matrices in object
#' 
#' The main input is the \code{tempoSig} object with count matrix.
#' This function performs non-negative factorization using Bayesian algorithm
#' and gamma priors. Slots \code{basis}, \code{coeff}, and \code{ranks} 
#' are filled.
#' 
#' @param object \code{scNMFSet} object containing count matrix.
#' @param ranks Rank for factorization; can be a vector of multiple values.
#' @param nrun No. of runs with different initial guesses.
#' @param progress.bar Display progress bar with \code{verbose = 1} for 
#'       multiple runs.
#' @param initializer If \code{'random'}, randomized initial conditions; 
#'        \code{'svd2'} for singular value decomposed initial condition.
#' @param Itmax Maximum no. of iteration.
#' @param hyper.update Vector of four logicals, each indcating whether
#'        hyperparameters \code{c(aw, bw, ah, bh)} should be optimized.
#' @param gamma.a Gamma distribution shape parameter.
#' @param gamma.b Gamma distribution mean. These two parameters are used for 
#'        fixed hyperparameters with \code{hyper.update} elements \code{FALSE}.
#' @param Tol Tolerance for terminating iteration.
#' @param hyper.update.n0 Initial number of steps in which hyperparameters 
#'        are fixed.
#' @param hyper.update.dn Step intervals for hyperparameter updates.
#' @param fudge Small positive number used as lower bound for factor matrix 
#'        elements to avoid singularity. If \code{fudge = NULL} (default), 
#'        it will be replaced by \code{.Machine$double.eps}. 
#'        Can be set to 0 to skip 
#'        regularization.
#' @param unif.stop Terminate if any of columns in basis matrix is uniform.
#' @return Object of class \code{scNMFSet} with factorization slots filled.
#' 
#' @details When run with multiple values of \code{ranks}, factorization is 
#'        repeated for each rank and the slot \code{measure} contains 
#'        log evidence and optimal hyperparameters for each rank. 
#'        With \code{nrun > 1}, the solution
#'        with the maximum log evidence is stored for a given rank.
#' @import Rcpp
#' @export
bnmf <- function(object, ranks=2:10, nrun=1, verbose=2, 
                         progress.bar=TRUE, initializer='random',
                         Itmax=10000, hyper.update=rep(TRUE,4), 
                         hyper.tying='ikj', normalize.sig=TRUE,
                         gamma.a=1, gamma.b=1, Tol=1e-5, 
                         hyper.update.n0=10, hyper.update.dn=1, 
                         fudge=NULL, kstar = 'kmax', useC = FALSE,
                         unif.stop=TRUE, sindex=NULL){
  
   if(!kstar %in% c('kmax','kopt')) stop('Unknown option for kstar')
   if(is.null(fudge)) fudge <- .Machine$double.eps
   mat <- catalog(object) # S4 class scNMFSet
   
   nullr <- sum(Matrix::rowSums(mat)==0)
   nullc <- sum(Matrix::colSums(mat)==0)
   if(nullr>0) stop('Input matrix contains empty rows')
   if(nullc>0) stop('Input matrix contains empty columns')
   
   ranks <- ranks[ranks <= ncol(mat)] # rank <= no. of columns
   nrank <- length(ranks)

   bundle <- list(mat=mat, ranks=ranks, verbose=verbose, gamma.a=gamma.a,
                  gamma.b=gamma.b, initializer=initializer, 
                  Itmax=Itmax, fudge=fudge, 
                  hyper.tying=hyper.tying,
                  hyper.update=hyper.update, 
                  hyper.update.n0=hyper.update.n0, 
                  hyper.update.dn=hyper.update.dn, Tol=Tol,
                  unif.stop=unif.stop, nrun=nrun, useC=useC, 
                  w=signat(object), h=expos(object), sindex=sindex)
   vb <- lapply(seq_len(nrun), FUN=bnmf_iterate, bundle)
  
   basis <- coeff <- vector('list',nrank)
   rdat <- c()
   ranks2 <- c()
   
   for(k in seq_len(nrank)){     # find maximum solutions for each rank
     rmax <- -Inf
     for(i in seq_len(nrun)){
       if(vb[[i]]$rdat[k] > rmax){
         imax <- i
         rmax <- vb[[i]]$rdat[k]
       }
     }
     if(rmax==-Inf) next
     ranks2 <- c(ranks2,ranks[k])
     rdat <- c(rdat,rmax)
     w <- vb[[imax]]$wdat[[k]]
     h <- vb[[imax]]$hdat[[k]]
     dw <- vb[[imax]]$dwdat[[k]]
     dh <- vb[[imax]]$dhdat[[k]]
     if(normalize.sig){
        cw <- colSums(w)
        w <- t(t(w) / cw)
        h <- h * cw
        dw <- t(t(dw) / cw^2)
        dh <- h * cw^2
     }
     basis[[k]] <- w
     coeff[[k]] <- h

#     if(initializer!='restart'){
       rownames(basis[[k]]) <- rownames(mat)
       colnames(coeff[[k]]) <- colnames(mat)
#     }
   }
   
   if(kstar=='kmax') Kstar <- max(which(is.finite(rdat)))
   else if(kstar=='kopt') Kstar <- which.max(rdat)
   
   signat(object) <- basis[[Kstar]]
   H <- coeff[[Kstar]]
   expos(object) <- t(t(H)/colSums(H))
   misc(object)[['measure']] <- data.frame(rank=ranks2, lml=rdat)

   return(object)
}

bnmf_iterate <- function(irun, bundle){
  
   nrow <- dim(bundle$mat)[1]
   ncol <- dim(bundle$mat)[2]
   nrank <- length(bundle$ranks)
   
   rdat <- rep(-Inf, nrank)
   wdat <- hdat <- dwdat <- dhdat <- hyperp <- list()

   if(bundle$verbose >= 2) if(bundle$nrun > 1) 
     cat('Run ',irun,'\n',sep='')
   
   for(irank in seq_len(nrank)){
     
     rank <- bundle$ranks[[irank]]
     if(rank > min(nrow,ncol))
       stop('K exceeded min(nrow,ncol)')
     aw <- bundle$gamma.a[1]
     ah <- bundle$gamma.a[length(bundle$gamma.a)]
     bw <- bundle$gamma.b[1]
     bh <- bundle$gamma.b[length(bundle$gamma.b)]
     
     hyper <- hyper0 <- list(aw=aw, bw=bw, ah=ah, bh=bh)
     if(bundle$initializer=='random')
       wh <- vb_init(nrow, ncol, bundle$mat, rank, hyper=hyper, 
                   initializer=bundle$initializer)
     else{
       if(!is.null(bundle$sindex))
         bundle$h[,bundle$sindex] <- stats::rgamma(n=rank, shape=ah)
       dw <- matrix(0, nrow=nrow, ncol=rank)
       dh <- matrix(0, nrow=rank, ncol=ncol)
       rownames(dw) <- rownames(bundle$mat)
       colnames(dw) <- seq_len(rank)
       rownames(dh) <- seq_len(rank)
       colnames(dh) <- colnames(bundle$mat)
       wh <- list(w=bundle$w, h=bundle$h, lw=bundle$w, lh=bundle$h, ew=bundle$w, 
                  eh=bundle$h, dw=dw, dh=dh)
     }
     lk0 <- 0
     for(it in seq_len(bundle$Itmax)){
       if(bundle$useC)
          wh <- vbnmf_update(as.matrix(bundle$mat), wh, hyper, c(bundle$fudge))
       else
          wh <- vbnmf_updateR(bundle$mat, wh, rank, hyper, fudge=bundle$fudge)
       if(bundle$initializer=='restart' | 
          it > bundle$hyper.update.n0 & it%%bundle$hyper.update.dn==0) 
         hyper <- hyper_update(bundle$hyper.update, wh, hyper, 
                               hyper.tying=bundle$hyper.tying, Niter=100, Tol=1e-3)
       if(is.na(wh$lkh)) break
       if(it>1) if(it > bundle$hyper.update.n0)
          if(wh$lkh>=lk0) if(abs(1-wh$lkh/lk0) < bundle$Tol) break
       lk0 <- wh$lkh
       if(bundle$verbose >= 3) cat(it,', lkl = ',lk0, '\n',sep='')
     }
     contains.unif <- apply(wh$ew,2,
                              function(x){abs(max(x)-min(x))<bundle$Tol})
     if(sum(contains.unif)>0){
#       cat('K ',rank,' row/column ',
#              paste(which(contains.unif),collapse=','),' constant.\n')
       cat('Kmax = ',rank - 1, '\n')
       break
     }
     rdat[irank] <- lk0
     wdat[[irank]] <- wh$ew
     hdat[[irank]] <- wh$eh
     dwdat[[irank]] <- sqrt(wh$dw)
     dhdat[[irank]] <- sqrt(wh$dh) 
     hyperp[[irank]] <- hyper
     if(bundle$verbose >= 2)
       cat('K = ',rank,': iteration = ',it,', lkl= ',lk0,'\n',sep='')
   }   # end of irank-loop
   
   vb <- list(rdat=rdat, wdat=wdat, hdat=hdat, hyperp=hyperp, dwdat=dwdat,
              dhdat=dhdat)
   return(vb)
}

