#' Determine optimal rank
#' 
#' @param object \code{scNMFSet} object containing factorization output, or 
#'        data frame containing the rank-evidence profile.
#' @param df Degrees of freedom for split fit. Upper bound is the total number of
#'        data points (number of rank values scanned).
#' @param BF.threshold Bayes factor threshold for statistical threshold. 
#' @param type \code{c(1,2)}. Type 1 is where there is a clear maximum. Type 2 
#'        is where marginal likelihood reaches a maximal level and stays constant.
#'        If omitted, the type will be inferred from data.
#' @param m Number of features (e.g., genes) in the count matrix. Only necessary when
#'        \code{object} is of type \code{data.frame}.
#' @return List containing \code{type} and \code{ropt} (optimal rank).
#' 
#' @details The input object is used along with Bayes factor threshold to determine the
#'          heterogeneity type (1 or 2) and the optimal rank. 
#'          If evidence(rank 1)/evidence(rank2) > \code{BF.treshold}, rank 1 is favorable than rank 2.
#' @examples
#' set.seed(1)
#' @export
kstar <- function(object, df=NULL, BF.threshold=3, type=1, m = 96){
  
  if(is(object,'tempoSig')){
    me <- misc(object)$measure[,c(1,2)]
    m <- nrow(catalog(object))
  }
  else if(is(object,'data.frame')){
    me <- object[,c(1,2)]
    if(is.null(m)) stop('No. of rows unknown')
  }
  else stop('Inappropriate class of object')
  
  if(is.null(df)) df <- nrow(me)
  df <- min(df,nrow(me))
  fs <- stats::smooth.spline(x=me[,1],y=me[,2],df=df)
  rst <- fs$x[which.max(fs$y)]     # arg max_r (L)
  bf <- log(BF.threshold)/m
  
#  if(is.null(type)){    # determine the type
#    rmax <- max(me[,1])
#    flag <- abs(fs$y-max(fs$y)) <= bf
#    range <- fs$x[flag]  # rank values within the threshold around optimum
#    if(rmax %in% range) type <- 2
#    else type <- 1
#  }
  
  if(type==1) ropt <- rst 
  else{
    sl <- slope(fs$y,fs$x)
    if(sum(sl < bf)>0) 
      idx <- which(sl < bf)[1]
    else
      idx <- nrow(me)
    ropt <- fs$x[idx]
  }
  
  return(list(type=type, ropt=ropt, spline=fs))
}

slope <- function(y,x){
  
  n <- length(x)
  s <- rep(0,n)
  s[1] <- (y[2]-y[1])/(x[2]-x[1])
  for(i in seq(2,n-1))
    s[i] <- (y[i+1]-y[i])/(x[i+1]-x[i])
  s[n] <- s[n-1]
  
  return(s)
}

#' Plot signature profile
#' 
#' @param x Signature profile vector; must be of length 96 and named for nucleotide contexts
#' @param lwd Line thickness for error bars
#' @return Return vector of \code{barplot}
#' @export
sigplot <- function(x, error = NULL, lwd = 0.5){
  
  nt <- trinucleotides()
  if(sum(is.na(names(x))) > 0 | sum(!names(x) %in% nt) > 0) 
    stop('Names of x not in trinucleotides')
  
  col0 <- c('cyan','black','red','tan','limegreen','pink')
  col <- rep(col0,each=16)
  
  x <- as.matrix(x)
  nt <- trinucleotides()
  if(is.null(error)) ylim <- c(0, max(c(x)))
  else ylim <- c(0, max(c(x) + c(error)))
  xt <- barplot(c(x),col=col,names.arg=rep('',96),xaxt='n',las=1,lwd=0.5, ylim=ylim)
  if(!is.null(error)){
    dx <- xt[2] - xt[1]
    segments(x0=xt, x1=xt, y0=c(x), y1=c(x) + c(error), col=col, lwd=lwd)
    segments(x0=xt-dx/3, x1=xt+dx/3, y0=c(x) + c(error), y1=c(x) + c(error), col=col, lwd=lwd)
  }
  ymax <- max(x)
  rect(xleft=xt[-96],xright=xt[-1],ytop=0,ybottom=-ymax*0.02,col=col,xpd=NA,lwd=0)
  axis(side=1,at=xt[seq(8,88,length.out=6)],label=c('C>A','C>G','C>T','T>A','T>C','T>G'),
       mgp=c(2,0.3,0),tck=0,lwd=0)
  title(main=colnames(x),cex.main=1.0)
  return(invisible(xt))
}

#' Filter exposures using p-values
#' 
#' @param object Object of class \code{tempoSig}
#' @param alpha Type-I error treshold
#' @param attribution If \code{TRUE}, output is multiplied by mutation counts of each sample
#' @return Matrix of filtered exposure or attribution
#' @export
filterExposure <- function(object, alpha = 0.05, attribution = TRUE){
  
  if(!is(object, 'tempoSig')) stop('Object is not of class tempoSig')
  E <- expos(object)
  pv <- pvalue(object)
  if(NROW(pv) == 0 | NCOL(pv) == 0) stop('p-value has not been filled')
  
  E2 <- E
  for(k in seq(NROW(E))) for(i in seq(NCOL(E)))
    if(pv[k, i] >= alpha) E2[k, i] <- 0
  if(attribution) E2 <- E2 * tmb(object)
  
  return(E2)
}
