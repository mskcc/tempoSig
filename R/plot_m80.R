#' Plot M80 Estimate
#' 
#' @param m80 Data frame outpuf from \code{powerEstimate}
#' @export
plot_M80 <- function(m80){
  
  sbs <- colnames(read.table(system.file('extdata', 'cosmic_SigAnalyzer_SBS_signatures.txt', 
                                package = 'tempoSig'), header = TRUE, sep = '\t'))
  Mmax <- ceiling(log10(max(m80$M80)))
  Signature <- paste(rownames(m80), m80$Etiology, sep=': ')
  dat1 <- data.frame(Signature = Signature, M80 = log10(m80$M80))
  dat2 <- data.frame(Signature = Signature, M80 = Mmax - log10(m80$M8))
  dat <- cbind(data.frame(box = c(rep('L',NROW(dat1)), rep('H',NROW(dat2)))), rbind(dat1,dat2))
  idx <- match(sapply(dat[,2], function(x){strsplit(as.character(x), split = ':')[[1]][1]}), sbs)
  expr <- c(1, 10, expression(10^2), expression(10^3), expression(10^4), expression(10^5),
            expression(10^6), expression(10^7))
  
  p <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x = M80, 
          y = reorder(Signature, NROW(dat1) - idx))) + 
    ggplot2::geom_bar(stat = 'identity', col = 'gray40', 
                      fill = c(rep(NA, NROW(dat1)), rep('gray60', NROW(dat2)))) +
    ggplot2::scale_y_discrete(position = 'right') + 
#   ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7, hjust = 1.0, vjust = 0.5)) + 
    ggplot2::scale_x_continuous(limits = c(0, Mmax), breaks = seq(0, Mmax), 
                                labels = expr[seq(Mmax + 1)]) + ggplot2::ylab('')
  print(p)
}
