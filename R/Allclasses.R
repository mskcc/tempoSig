#' Class \code{tempoSig} for storing catalog and exposure
#'
#' \code{S4} class with storage for input catalog, reference signatures,
#' and output exposure inferred.
#'
#' @slot catalog Matrix of mutation count with samples in columns
#'               and mutation types in rows
#' @slot signat Matrix of reference signatures; each column
#'                 stores individual signatures proportions with mutation
#'                 types in rows
#' @slot tmb Vector of total mutation loads for each sample
#' @slot expos Matrix of inferred proportions of signatures
#'                exposure; each row gives proportions of each
#'                sample with signatures in columns
#' @slot pvalue P-values of exposures estimated from permutation tests;
#'            same dimension as \code{expos}.
#' @return Object of class \code{tempoSig}
#' @export tempoSig
setClass('tempoSig',
         slots = c(catalog = 'matrix',
                   signat = 'matrix',
                   tmb = 'vector',
                   expos = 'matrix',
                   pvalue = 'matrix'
         ))

#' Create \code{tempoSig} object
#'
#' @param data Input data of mutation catalog; samples in columns
#'             and mutation types in rows.
#' @param signatures Reference signatures table; default is COSMIC
#'                  version 3 signatures.
#' @return Object of class \code{tempoSig}.
#' @examples
#' set.seed(130)
#' data <- rmultinom(n = 5, size = 100, prob = rep(1,96))
#' rownames(data) <- trinucleotides()
#' colnames(data) <- paste0('S',seq(5))
#' s <- tempoSig(data)
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt',
#'                   package='tempoSig'), header=TRUE, sep='\t')
#' b <- tempoSig(data)
#' @export
tempoSig <- function(data, signat = NULL){

  if(!is(data, 'matrix')) data <- as.matrix(data)
  if(is.null(signat))
    signat <- as.matrix(read.table(system.file('extdata/cosmic_sigProfiler_SBS_signatures.txt',
                            package = 'tempoSig')))
  else if(is.character(signat)){
    if(!file.exists(signat)) stop(paste0('File ',signat,' does not exist'))
    signat <- as.matrix(read.table(signat))
  } else if(!is(signat, 'matrix')) signat <- as.matrix(signat)
  nts <- rownames(signat)
  ntd <- rownames(data)
  if(sum(is.na(ntd)) > 0) stop('Data must have explicit mutation type names')
  idx <- match(nts, ntd)
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Row names in data do not match reference signatures')
  data <- data[idx, , drop = FALSE]

  x <- new('tempoSig', catalog = data, signat = signat)
  x@tmb <- colSums(data)

  return(x)
}

#' Display \code{tempoSig} object
#'
#' Display the class and dimension
#'
#' @param object Object of class \code{tempoSig}
#' @return \code{NULL}
#' @export
setMethod('show', signature = 'tempoSig',
          definition = function(object){
            cat('Class:', class(object), '\n')
            ca <- catalog(object)
            cat('Number of mutation types = ', NROW(ca), '\n')
            cat('  ')
            print(utils::head(rownames(ca)))
            cat('Number of samples = ', NCOL(ca), '\n')
            sg <- signat(object)
            cat('Number of signatures = ', NCOL(sg), '\n')
            cat('  ')
            print(utils::head(colnames(sg)))
          })

#' @export
setGeneric('catalog', function(object) standardGeneric('catalog'))
#' Accessor for catalog
#'
#' @param object Object containing catalog
#' @return \code{catalog} data frame
#' @export
setMethod('catalog', signature = 'tempoSig',
          function(object){
            object@catalog
          }
)
#' @export
setGeneric('signat', function(object) standardGeneric('signat'))
#' Accessor for signature
#'
#' @param object Object containing signatures
#' @return \code{signat} data frame
#' @export
setMethod('signat', signature = 'tempoSig',
          function(object){
            object@signat
          }
)
#' @export
setGeneric('tmb', function(object) standardGeneric('tmb'))
#' Accessor for TMB
#'
#' @param object Object containing \code{tmb}
#' @return \code{tmb} data frame
#' @export
setMethod('tmb', signature = 'tempoSig',
          function(object){
            object@tmb
          }
)
#' @export
setGeneric('tmb<-', function(object, value) standardGeneric('tmb<-'))
#' @export
setMethod('tmb<-', signature = 'tempoSig',
          function(object, value){
            object@tmb <- value
            if(validObject(object)) return(object)
          })
#' @export
setGeneric('expos', function(object) standardGeneric('expos'))
#' Accessor for exposure
#'
#' @param object Object containing \code{expo}
#' @return \code{expos} data frame
#' @export
setMethod('expos', signature = 'tempoSig',
          function(object){
            object@expos
          }
)
#' @export
setGeneric('expos<-', function(object, value) standardGeneric('expos<-'))
#' @export
setMethod('expos<-', signature = 'tempoSig',
          function(object, value){
            object@expos <- value
            if(validObject(object)) return(object)
})
#' @export
setGeneric('pvalue', function(object) standardGeneric('pvalue'))
#' Accessor for p-values
#'
#' @param object Object containing \code{pvalue}
#' @return \code{pvalue} data frame
#' @export
setMethod('pvalue', signature = 'tempoSig',
          function(object){
            object@pvalue
          }
)
#' @export
setGeneric('pvalue<-', function(object, value) standardGeneric('pvalue<-'))
#' @export
setMethod('pvalue<-', signature = 'tempoSig',
          function(object, value){
            object@pvalue <- value
            if(validObject(object)) return(object)
          })
