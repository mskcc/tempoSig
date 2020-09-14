#' Wrapper Function
#' 
#' One-step execution of exposure inference, filtering, and output
#' 
#' @param catalog Catalog input file
#' @param signat Reference signature choices, \code{c('v2','SA','SP')} for sCOSMIC v2, 
#'        SigAnalyzer, and sigProfiler, respectively
#' @param nperm No. of permutations 
#' @param alpha False positive rate threshold
#' @param output Output file
#' @export
run_tSig <- function(catalog, signat = 'v2', nperm = 1000, alpha = 0.05, output = 'exposure.txt'){
  
  if(signat=='SA') ref <- system.file('extdata', 'cosmic_SigAnalyzer_SBS_signatures.txt', package = 'tempoSig')
  else if(signat=='SP') ref <- system.file('extdata', 'cosmic_sigProfiler_SBS_signatures.txt', package = 'tempoSig')
  else if(signat=='v2') ref <- system.file('extdata', 'cosmic_snv_signatures_v2.txt', package = 'tempoSig')
  else stop('Unknown choice for signat')
  b <- tempoSig(catalog, signat = ref)
  b <- extractSig(b, compute.pval = TRUE, nperm = nperm, progress.bar = TRUE)
  e <- expos(b)
  pv <- pvalue(b)
  for(i in seq(NROW(e))) for(j in seq(1, NCOL(e))){
    if(is.na(e[i,j])) next()
    if(pv[i,j] >= alpha) e[i,j] <- 0
  }
  e <- cbind(data.frame(Sample.Name = rownames(e), Number.of.Mutations = tmb(b)), e)
  write.table(e, file = output, row.names = F, sep = '\t', quote=F)
}