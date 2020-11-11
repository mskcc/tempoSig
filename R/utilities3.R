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
kstar <- function(object, df=10, BF.threshold=3, type=1, m = 96){
  
  if(is(object,'tempSig')){
    me <- misc(object)$measure[,c(1,2)]
    m <- nrow(catalog(object))
  }
  else if(is(object,'data.frame')){
    me <- object[,c(1,2)]
    if(is.null(m)) stop('No. of rows unknown')
  }
  else stop('Inappropriate class of object')
  
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
  
  return(list(type=type, ropt=ropt))
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
#' @param ... Other parameters for \code{plot}
#' @return Return vector of \code{barplot}
#' @export
sigplot <- function(x, ...){
  
  nt <- trinucleotides()
  if(sum(is.na(names(x))) > 0 | sum(!names(x) %in% nt) > 0) 
    stop('Names of x not in trinucleotides')
  
  col0 <- c('cyan','black','red','tan','limegreen','pink')
  col <- rep(col0,each=16)
  
  x <- as.matrix(x)
  nt <- trinucleotides()
  xt <- barplot(c(x),col=col,names.arg=rep('',96),xaxt='n',las=1,lwd=0.5)
  ymax <- max(x)
  rect(xleft=xt[-96],xright=xt[-1],ytop=0,ybottom=-ymax*0.02,col=col,xpd=NA,lwd=0)
  axis(side=1,at=xt[seq(8,88,length.out=6)],label=c('C>A','C>G','C>T','T>A','T>C','T>G'),
       mgp=c(2,0.3,0),tck=0,lwd=0)
  title(main=colnames(x),cex.main=1.0)
  return(xt)
}

