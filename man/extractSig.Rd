% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.R
\name{extractSig}
\alias{extractSig}
\title{Infer Signature Proportions}
\usage{
extractSig(
  object,
  method = "mle",
  itmax = 1000,
  tol = 1e-04,
  min.tmb = 2,
  compute.pval = FALSE,
  nperm = 1000,
  progress.bar = FALSE,
  pvtest = "permutation",
  cosmic = FALSE,
  Kmin = 2,
  K = 2,
  Kmax = 30,
  nrun = 10,
  useC = TRUE,
  initializer = "random",
  verbose = 2,
  ...
)
}
\arguments{
\item{object}{Object of class \code{tempoSig}}

\item{method}{Refitting method; \code{mle} for maximum likelihood (default) or
\code{mutCone} for mutationalCone.}

\item{itmax}{Maximum number of iterations for maximum likelihood estimate}

\item{tol}{Tolerance for convergence}

\item{min.tmb}{Minimum number of mutations in each sample. If \code{tmb} is less,
\code{NA} will be returned for the sample.}

\item{compute.pval}{Estimate p-values}

\item{nperm}{Number of permutations}

\item{progress.bar}{Display progress bar}

\item{pvtest}{Algorithm for p-value computation; 
\code{c('permutation','lrt','x.permutation')} for permutation resampling of 
signatures, likelihood ratio test (asymptotic formula), or
permutation of count data.}

\item{...}{Other parameters for \code{denovo} with \code{method = 'hnmf'}}
}
\description{
Use known signature list to find the most likely exposures in samples
}
\examples{
data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt', package='tempoSig'))
b <- tempoSig(data)
b <- extractSig(b, progress.bar = TRUE)
b_pv <- extractSig(b, compute.pval = TRUE, progress.bar = TRUE)
}
