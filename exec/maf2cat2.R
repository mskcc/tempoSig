#!/usr/bin/env Rscript
# Convert maf file with Ref_Tri column into catalog
library('argparse')
library('tempoSig', lob.loc = '/usr/local/lib/R/site-library')

parser <- ArgumentParser(description = 'Construct mutational catalog from MAF file with Ref_Tri column')
parser$add_argument('maf', metavar = 'MAF',
                    help = 'input MAF file')
parser$add_argument('catalog', metavar = 'CATALOG', 
                    help = 'output catalog file')
args <- parser$parse_args()

catalog <- maf2cat(maf = args$maf)
write.table(catalog, file = args$catalog, quote = F, sep = '\t')
