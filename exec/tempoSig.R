#!/usr/bin/env Rscript
library('argparse')
library('tempoSig', lib.loc = '/usr/local/lib/R/site-library')
parser <- ArgumentParser(description = 'Fit mutational catalog to signatures')
parser$add_argument('catalog', metavar = 'CATALOG', 
                    help = 'input catalog data file')
parser$add_argument('output', metavar = 'OUTPUT', help = 'output file name')
group <- parser$add_mutually_exclusive_group()
group$add_argument('--cosmic_v2', action = 'store_true',
                   help = 'use COSMIC v2 reference signatures (default)')
group$add_argument('--cosmic_v3', action = 'store_true', 
                   help = 'use COSMIC v3 reference signatures')
group$add_argument('--cosmic_v3_SA', action = 'store_true',
                   help = 'use COSMIC v3 SigAnalyzer reference signatures')
group$add_argument('--cosmic_v3_exome', action = 'store_true',
                   help = 'use COSMIC v3 exome reference signatures')
parser$add_argument('--sigfile', dest = 'sigfile', action = 'store',
                    help = 'custom input reference signature file; overrides --cosmic_v2/3')
parser$add_argument('--pvalue', action = 'store_true', default = FALSE,
                    help = 'estimate p-values (default FALSE)')
parser$add_argument('--nperm', dest = 'nperm', action = 'store', type = 'integer',
                    default = 1000, help = 'number of permutations for p-value estimation; default 1000')
parser$add_argument('--seed', dest = 'seed', action = 'store', type = 'integer',
                    help = 'random number seed')
parser$add_argument('--pv.out', dest = 'pv.out', action = 'store', 
                    help = 'p-value output file')
parser$add_argument('--cbio', action = 'store_true', default = FALSE,
                    help = 'output in cBioPortal format (default FALSE)')
args <- parser$parse_args()

if(!is.null(args$seed)) set.seed(args$seed)

if(!is.null(args$sigfile)){  # custom signature file provided
  if(!file.exists(args$sigfile)) stop(paste0(args$sigfile, ' does not exist.'))
  sig <- args$sigfile                                      
} else if(args$cosmic_v3){
  sig <- 'SP'
} else if(args$cosmic_v3_exome){
  sig <- system.file('extdata', 'cosmic_sigProfiler_exome_SBS_signatures.txt', package = 'tempoSig')
} else if(args$cosmic_v3_SA){
  sig <- 'SA'
} else{  # cosmic_v2 (default)
  sig <- 'v2'
}

if(args$cbio & args$pvalue & is.null(args$pv.out)){
  args$pv.out <- 'pvalues.txt'
}

data <- read.table(args$catalog, header = TRUE, sep = '\t')
x <- tempoSig(data = data, signat = sig)
x <- extractSig(object = x, compute.pval = args$pvalue, nperm = args$nperm, progress.bar = TRUE)
writeExposure(object = x, output = args$output, pv.out = args$pv.out, sep = '\t', 
              cBio.format = args$cbio)
