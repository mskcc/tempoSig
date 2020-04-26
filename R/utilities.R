#' Trinucleotide labels
#'
#' @export
trinucleotides <- function(nt = 'ACGT', brackets = TRUE, arrows = TRUE){

  nt2 <- strsplit(nt, split='')[[1]]
  pyr <- nt2[nt2 %in% c('C', 'T')]
  w <- expand.grid(nt2, pyr, stringsAsFactors = FALSE)
  w <- w[,c(2,1)]
  w <- w[w[, 1]!=w[,2],]
  sp <- ifelse(arrows, '>', '')
  wx <- apply(w, 1, function(x){paste(c(x[1], sp, x[2]),collapse='')})

  z <- expand.grid(nt2, nt2, wx, stringsAsFactors = FALSE)
  z <- z[, c(2, 3, 1)]
  if(brackets) sp <- c('[',']')
  else sp <- c('','')
  w <- apply(z, 1, function(x){paste(c(x[1], sp[1], x[2], sp[2], x[3]), collapse = '')})
  names(w) <- NULL
  return(w)
}

#' Plot exposure
#'
#' Display barplot of exposure extracted for a specific sample
#'
#' @param object Object of class \code{tempoSig}
#' @param sample.id Sample index or name (\code{Tumor_Sample_Barcode}) to display
#' @param cutoff Minimum proportion for displaying signature labels
#' @param ... Other parameters to \code{barplot}
#' @return None
#' @export
plotExposure <- function(object, sample.id, cutoff = 1e-3, ...){

  if(!is(object, 'tempoSig')) stop('Object is not of class tempoSig')
  sname <- names(tmb(object))
  if(is.character(sample.id)){
    sample.id <- which(sname == sample.id)
    if(length(sample.id)==0) stop(paste0(sample.id, 'is not in object'))
  }
  sample.id <- as.integer(sample.id)
  expo <- exposure(object)
  if(sample.id < 1 | sample.id > NROW(expo)) stop('sample.id out of bound in object')
  e <- expo[sample.id,]
  names.arg <- names(e)
  names.arg[e < cutoff] <- ''
  graphics::barplot(e, main = sname[sample.id], las=2, names.arg = names.arg,
                    ylab = 'Proportions', ...)

  return(invisible(e))
}

#' Write Exposure
#'
#' Save a text output of exposures
#'
#' Writes a text file of specified name.
#'
#' @param object Object of class \code{tempoSig}
#' @param output File name of the output
#' @param sep Delimiter, either space or tab.
#'
#' @export
writeExposure <- function(object, output, sep = '\t'){

  if(!is(object, 'tempoSig')) stop('Object is not of class tempoSig')
  if(!is.character(output)) stop('Output file name must be characters')
  if(!sep %in% c(' ','\t')) stop('Delimiter must be either space or tab')
  expo <- exposure(object)
  if(all(dim(expo) == 0)) stop('Exposure in object empty')

  is.pv <- !all(dim(pvalue(object)) == 0)   # pvalue is not empty
  if(!is.pv)
    out <- cbind(data.frame(Tumor_Sample_Barcode = rownames(expo), TMB = tmb(object)),
               as.data.frame(expo))
  else{
    out <- data.frame(Tumor_Sample_Barcode = rownames(expo), TMB = tmb(object))
    pv <- pvalue(object)
    sig.names <- colnames(expo)
    for(k in seq(NCOL(expo))){
      tmp <- data.frame(expo[,k], pv[,k])
      names(tmp) <- paste0(sig.names[k], c('.observed','.pvalue'))
      out <- cbind(out, tmp)
    }
  }
  write.table(out, file=output, sep = sep, row.names = F, quote = F)
  return(invisible(object))
}