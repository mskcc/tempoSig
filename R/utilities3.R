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

#' Construct catalog matix from MAF file (without of Ref_Tri column)
#' 
#' @param maf MAF file or data frame of chromosome, position, ref, and alt alleles
#' @param ref.genome Reference genome of class \code{BSgenome}
#' @param fix.chr Modify chromosome names into \code{chr1, chr2, ...} for \code{ref.genome}
#' @param progress.bar Display progress bar
#' @return Catalog matrix
#' @examples 
#' library(BSGenome.Hsapiens.UCSC.hg19)
#' maf <- system.file('extdata', 'brca.maf', package = 'tempoSig')
#' x <- maf2cat3(maf, BSGenome.Hsapiens.HCSC.hg19)
#' head(x)
#' 
#' @export
maf2cat3 <- function(maf, ref.genome, fix.chr = TRUE, progress.bar = TRUE){
  
  if(is(maf, 'character')){
    if(!file.exists(maf)) stop(paste0('File ', maf, ' does not exist'))
    maf <- data.table::fread(maf, data.table = F)
  } 
  if(fix.chr) maf <- readMAFforcat(maf = maf)
  
  colnames(maf)[colnames(maf) == 'Start_position'] <- 'Start_Position'
  colnames(maf)[colnames(maf) == 'End_position'] <- 'End_Position'
  
  if('Variant_Type' %in% colnames(maf)) 
    maf <- maf[maf$Variant_Type == 'SNP', ]
  
  cmplmnt <- c('A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')
  
  coln <- c('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode')
  if(sum(!coln %in% colnames(maf)) > 0) 
    stop(paste0('Necessary columns of MAF are missing: ', coln, collapse = ''))

  if(!is(ref.genome, 'BSgenome'))
     stop('ref.genome must be of class BSgenome')
  
  chr <- unique(maf[,'Chromosome'])
  if(sum(!chr %in% seqnames(ref.genome)) > 0) stop('There is an unknown chromosome name')
  
  ref.chr <- seqnames(ref.genome)
  nchr <- length(chr)
  
  if(progress.bar)
    pb <- txtProgressBar(style = 3)

  nt <- trinucleotides()
  tsb <- unique(maf$Tumor_Sample_Barcode)
  nsample <- length(tsb)
  xcat <- matrix(0, nrow = length(nt), ncol = nsample)
  rownames(xcat) <- nt
  colnames(xcat) <- tsb
  
  for(k in seq_along(ref.chr)){
    kchr <- ref.chr[k]
    dat <- maf[maf$Chromosome==kchr, ]
    if(NROW(dat) == 0) next()
    rchr <- ref.genome[[kchr]]
    dat <- dat[dat$Start_Position > 1 & dat$Start_Position < length(rchr), ]
    m <- NROW(dat)
    if(m ==0) next()
    # subset of variants within reference positions
    pos <- dat$Start_Position
    ref.alleles <- rchr[pos]
    if(as.character(ref.alleles) != paste0(dat$Reference_Allele, collapse = ''))
      stop('Ref. alleles in data and ref.genome do not match')
    
    ref.tri <- matrix(NA, nrow = m, ncol = 4)
    ref.tri[, 2] <- dat$Reference_Allele
    ref.tri[, 3] <- dat$Tumor_Seq_Allele2
    ref.tri[, 1] <- strsplit(as.character(rchr[pos - 1]), split = '')[[1]]
    ref.tri[, 4] <- strsplit(as.character(rchr[pos + 1]), split = '')[[1]]
    
    purines <- ref.tri[, 2] %in% c('A', 'G')
    for(i in seq(4))
      ref.tri[purines, i] <- cmplmnt[ref.tri[purines, i]]
    ref.tri[purines, c(1, 4)] <- ref.tri[purines, c(4, 1)] # got to reverse the +-1 position
    ref.tri <- apply(ref.tri, 1, 
                     function(x){paste0(x[1], '[', x[2], '>', x[3], ']', x[4], collapse = '')})
    
    cnt <- as.matrix(table(ref.tri, dat$Tumor_Sample_Barcode))
    cnt <- cnt[rownames(cnt) %in% nt, , drop = F]  # discard mutation types not in 96 (eg., N[C>T]A)
    xcat[rownames(cnt), colnames(cnt)] <- xcat[rownames(cnt), colnames(cnt)] + cnt
    
    if(progress.bar) setTxtProgressBar(pb, k/nchr)
  }
  if(progress.bar) close(pb)
  
  return(xcat)
}

#' Read MAF file and Fix Chromosome Names
#' 
#' This function is a wrapper intended for use in \code{maf2cat3}, which requires chromosome names
#' compatible with \code{BSgenome.Hsapiens.UCSC.hg19} (\code{chr1, ..., chrX, chrY}).
#' 
#' @param maf MAF file path
#' @export
#' @return Data frame with chromosome names fixed
#' @examples
#' maf <- system.file('extdata', 'brca.maf', package = 'tempoSig')
#' maf <- readMAFforcat(maf)
#' head(maf)
#' 
readMAFforcat <- function(maf){
  
  if(is.character(maf)){
    if(!file.exists(maf)) stop(paste0('File ',maf, 'does not exist'))
    x <- data.table::fread(maf, data.table = F)
  } else
    x <- as.data.frame(maf)
  
  flag <- !grepl('chr', x$Chromosome)
  x[flag, 'Chromosome'] <- paste0('chr', x[flag, 'Chromosome'])
  x[x$Chromosome == 'chr23', 'Chromosome'] <- 'chrX'
  x[x$Chromosome == 'chr24', 'Chromosome'] <- 'chrY'
  
  return(x)
}