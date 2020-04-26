#!/usr/bin/env Rscript
library('argparse')
library(tempoSig)
parser <- ArgumentParser(description = 'Fit mutational catalog to signatures')
parser$add_argument('catalog', metavar = 'CATALOG', 
                    help = 'input catalog data file')
parser$add_argument('output', metavar = 'OUTPUT', help = 'output file name')
group <- parser$add_mutually_exclusive_group()
group$add_argument('--cosmic_v3', action = 'store_true', 
                   help = 'use COSMIC v3 reference signatures (default)')
group$add_argument('--cosmic_v2', action = 'store_true',
                   help = 'use COSMIC v2 reference signatures')
parser$add_argument('--sigfile', dest = 'sigfile', action = 'store',
                    help = 'custom input reference signature file; overrides --cosmic_v3/2')
parser$add_argument('--pvalue', action = 'store_true', default = FALSE,
                    help = 'estimate p-values')
parser$add_argument('--nperm', dest = 'nperm', action = 'store', type = 'integer',
                    default = 1000, help = 'number of permutations for p-value estimation; default 1000')
args <- parser$parse_args()

if(!is.null(args$sigfile)){  # custom signature file provided
  if(!file.exists(args$sigfile)) stop(paste0(args$sigfile, ' does not exist.'))
  sig <- args$sigfile                                      
} else if(args$cosmic_v2){
  sig <- system.file('extdata', 'COSMIC_SNV_signatures_v2.txt', package = 'tempoSig')
} else{  # cosmic_v3 (default)
  sig <- NULL
}

data <- read.table(args$catalog, header = TRUE, sep = '\t')
x <- tempoSig(data = data, signat = sig)
x <- extractSig(object = x, compute.pval = args$pvalue, nperm = args$nperm, progress.bar = TRUE)
writeExposure(object = x, output = args$output, sep = '\t')
