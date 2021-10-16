#!/usr/bin/env Rscript
library('argparse')
library('tempoSig')
parser <- ArgumentParser(description = 'Perform piggyback de novo inference')
parser$add_argument('catalog', metavar = 'CATALOG', 
                    help = 'input catalog data file')
parser$add_argument('output', metavar = 'OUTPUT', help = 'output file name')
group <- parser$add_mutually_exclusive_group()
group$add_argument('--filter', action = 'store_true', default = FALSE,
                   help = 'filter exposures using CV cutoffs (default FALSE)')
parser$add_argument('--cutoff', dest = 'cutoff', action = 'store',
                    help = 'custom cutoff file')
parser$add_argument('--seed', dest = 'seed', action = 'store', type = 'integer',
                    help = 'random number seed')
args <- parser$parse_args()

if(!is.null(args$seed)) set.seed(args$seed)

if(!is.null(args$cutoff)){  # custom signature file provided
  if(!file.exists(args$cutoff)) stop(paste0(args$cutoff, ' does not exist.'))
  cutoff <- args$cutoff
} else{  
  cutoff <- NULL
}

data <- read.table(args$catalog, header = TRUE, sep = '\t')
x <- pgnmf(x = data, progress.bar = TRUE, filter = args$filter, cutoff = args$cutoff)

write.table(x$mean, file=args$output, sep = '\t', quote = F, row.names = TRUE, col.names = TRUE)
