#' Mutation load prediction plot
#' @export
powerPlot <- function(mpredict){
  
  mpr <- mpredict[is.finite(mpredict[,2]), 2, drop = FALSE]
  mmax <- 10^ceiling(log10(max(mpr)))
  dat <- t(as.matrix(cbind(mpr, mmax - mpr)))
  mp <- barplot(dat, beside = FALSE, col = c('white', 'gray'), las=2, log='y', 
                ylab='Mutation load (power > 0.8)')
}