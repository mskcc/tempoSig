#' Infer Signature Proportions
#'
#' Use known signature list to find the most likely exposures in samples
#'
#' @param spectra Matrix of mutational spectra with mutation contexts in rows and genome
#'        samples in columns
#' @param ref Reference signature database with mutation contexts in rows and signatures
#'        in columns
#' @param maf Data frame of maf containing aggregated mutation data
#' @import Rcpp
#' @useDynLib tempoSig
#' @export
extractSig <- function(spectrum = NULL, maf = NULL, ref = NULL, itmax = 1000, tol = 1e-4){

  if(!is.null(spectrum))
    spectrum <- t(spectrum)  # input has genomes in rows and mutatino types in columns
  mut.load <- colSums(spectrum)

  if(is.null(spectrum) & is.null(maf))
    stop('Either spectrum or maf must be provided.')

  if(is.null(ref)){
#   fl <- system.file('extdata/cosmic_SBS.signatures_v3.txt', package = 'tempoSig')
    fl <- system.file('extdata/Stratton_signatures30.txt', package = 'tempoSig')
    ref <- read.table(fl, header = TRUE, sep = '\t')
    ref <- as.matrix(t(ref))  # rows: trinucleotides, columns: SBS
  }
  nt <- rownames(ref)
  idx <- match(nt, rownames(spectrum))
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Mutation types in spectrum do not match reference.')
  spectrum <- spectrum[idx, ]  # rearrange rows to match reference

  h <- apply(spectrum, 2, decompose, ref = ref, itmax = itmax, tol = tol)
  rownames(h) <- colnames(ref)

  H <- cbind(data.frame(Sample.Name = colnames(spectrum),
                  Number.of.Mutations = mut.load), t(h))
  return(H)
}

# x = vector of observed proportions of mutation contexts
# ref = matrix of reference signature proportions
decompose <- function(x, ref, itmax = 1000, tol = 1e-4){

  num_sigs <- NCOL(ref)
  num_muts <- sum(x)
  if(num_muts < 5){
    warning('Sample has less than 5 mutations.')
    return(rep(NA, num_muts))
  }
  x <- x / num_muts
  x0 <- gtools::rdirichlet(n = 1, alpha = rep(10, num_sigs)) #initial guess
  p <- mlestimate(x, x0, ref, Itmax=itmax, Tol=tol)
  h <- p$x^2/sum(p$x^2)

  return(h)
}
