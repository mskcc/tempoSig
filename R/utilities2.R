#' Filter False Positives
#' 
#' Use exposure and p-value tables and generated filtered exposures
#' 
#' @param exposure Exposure input data file
#' @param pval P-value input data file
#' @param alpha FPR threshold
#' @param min.M Minimum mutation count
#' @export
#' @return Data frame of filtered exposure data
#' 
filterFP <- function(exposure, pval, alpha = 0.05, min.M = NA){
  
  if(is(exposure,'character'))
    e <- read.table(exposure, header = TRUE, sep = '\t')
  else if(is(exposure,'data.frame'))
    e <- exposure
  else stop('Unknown input type for exposure')
  
  if(is(pval,'character'))
    p <- read.table(pval, header = TRUE, sep = '\t')
  else if(is(pval, 'data.frame'))
    p <- pval
  else stop('Unknown input type for pvalue')
  
  if(!all(colnames(e) == colnames(p))) stop('Colnames of exposure and pvalues do not match')
  if(!all(e[,'Sample.Name'] == p[,'Sample.Name'])) stop('Sample names of exposure and pvalues do not match')
  if(!is.na(min.M)) if(! ('Number.of.Mutations' %in% colnames(e))) stop('No. of mutations not in data')
  
  e2 <- e
  for(i in seq(nrow(e))){
    if(is.na(e[i,3])) next()
    if(!is.na(min.M) & e[i,'Number.of.Mutations'] < min.M){ 
      e2[i,seq(3,ncol(e))] <- NA
      next()
    }
    for(j in seq(3,ncol(e)))
      if(p[i,j] >= alpha) e2[i,j] <- 0
  }
  
  return(e2)
}